"""Tests for GroupFit: sequential fits over a group of scans."""

import numpy as np
import pytest

from tavi.library.data.scan import UUID, Provenance, RawScan, ScanData, ScanMetadata, TaviMetadata
from tavi.library.fit import FitResult, GroupFit, GroupFitResult, ModelName

S1 = np.linspace(-2.0, 3.0, 101)


def make_scan(
    uuid_val: str,
    temp: float,
    center: float,
    amplitude: float = 1000.0,
    sigma: float = 0.2,
    background: float = 10.0,
    monitor: float = 1.0,
) -> RawScan:
    """A Poisson-noisy Gaussian in detector vs s1 at a near-constant temperature, with a monitor column."""
    expected = background + amplitude / (sigma * np.sqrt(2 * np.pi)) * np.exp(-((S1 - center) ** 2) / (2 * sigma**2))
    # Real counts are noisy; noiseless data makes lmfit's covariance degenerate (no stderr).
    detector = np.random.default_rng(0).poisson(expected * monitor).astype(float)
    return RawScan(
        uuid=UUID(value=uuid_val),
        data=ScanData(
            data={
                "s1": S1.tolist(),
                "detector": detector.tolist(),
                "temp": (temp + 0.01 * np.sin(S1)).tolist(),
                "mcu": [monitor] * len(S1),
            }
        ),
        metadata=ScanMetadata(),
        tavimeta=TaviMetadata(default_axis=("s1", "detector"), friendly_name=uuid_val, friendly_path="/exp1"),
        prov=Provenance(raw_file=f"{uuid_val}.dat", contributing_scans={UUID(value=uuid_val): 1}),
    )


def make_sweep() -> list[RawScan]:
    """Three scans whose peak walks to the right as temperature rises."""
    return [make_scan("s-10", 10.0, 0.0), make_scan("s-20", 20.0, 0.3), make_scan("s-30", 30.0, 0.6)]


PEAK_AND_BACKGROUND = [
    (ModelName.Gaussian, dict(guess=True, prefix="peak1_")),
    (ModelName.Linear, dict(slope=0.0, intercept=10.0, prefix="bkg_")),
]


def test_fit_returns_one_result_per_scan_in_order() -> None:
    """Each scan gets a GroupFitResult carrying its own uuid, in the order given."""
    scans = make_sweep()
    results = GroupFit(scans).fit([(ModelName.Gaussian, dict(guess=True))])
    assert [r.uuid for r in results] == [s.uuid for s in scans]
    assert all(isinstance(r, GroupFitResult) and isinstance(r.fit_result, FitResult) for r in results)


def test_summarize_center_vs_temperature() -> None:
    """Fitted centers line up with the injected centers against the mean temperature column."""
    group = GroupFit(make_sweep())
    group.fit(PEAK_AND_BACKGROUND)
    summary = group.summarize("temp", "center")
    assert np.allclose(summary.x, [10.0, 20.0, 30.0], atol=0.01)
    assert np.all(summary.x_err > 0)
    assert np.allclose(summary.y, [0.0, 0.3, 0.6], atol=0.01)
    assert np.all(np.isfinite(summary.y_err))
    assert summary.uuids == [UUID(value="s-10"), UUID(value="s-20"), UUID(value="s-30")]


def test_summarize_by_prefix_and_chi_squared() -> None:
    """A prefix reads a non-peak component; reduced_chi_squared reads fit quality with NaN error."""
    group = GroupFit(make_sweep())
    group.fit(PEAK_AND_BACKGROUND)
    intercept = group.summarize("temp", "intercept", prefix="bkg_")
    assert np.allclose(intercept.y, 10.0, atol=2.0)
    chi = group.summarize("temp", "reduced_chi_squared")
    assert np.all(np.isfinite(chi.y)) and np.all(np.isnan(chi.y_err))


def test_norm_col_divides_intensity() -> None:
    """Normalizing by the monitor column undoes a monitor-scaled detector."""
    scans = [make_scan("a", 10.0, 0.0, monitor=4.0), make_scan("b", 20.0, 0.0, monitor=2.0)]
    group = GroupFit(scans)
    group.fit(PEAK_AND_BACKGROUND, norm_col="mcu")
    assert np.allclose(group.summarize("temp", "amplitude").y, 1000.0, rtol=0.05)


def test_x_range_masks_points() -> None:
    """Only points inside x_range are fit."""
    results = GroupFit(make_sweep()).fit([(ModelName.Gaussian, dict(guess=True))], x_range=(1.0, -1.0))
    n_points = len(results[0].fit_result.best_fit)
    assert n_points == np.count_nonzero((S1 >= -1.0) & (S1 <= 1.0))


def test_chain_seeds_from_previous_fit() -> None:
    """With chain, each scan starts from its neighbour's fit, so a narrow peak drifting away is tracked."""
    centers = [2.0 + 0.05 * i for i in range(8)]
    scans = [make_scan(f"s{i}", 10.0 * i, c, sigma=0.03) for i, c in enumerate(centers)]
    # The peak drifts ~12 widths from where every scan's fit starts unchained.
    start = [(ModelName.Gaussian, dict(center=2.0, sigma=0.03, amplitude=1000.0))]
    group = GroupFit(scans)
    group.fit(start)
    assert not np.allclose(group.summarize("temp", "center").y, centers, atol=0.01)
    group.fit(start, chain=True)
    assert np.allclose(group.summarize("temp", "center").y, centers, atol=0.01)


def test_chain_keeps_set_options() -> None:
    """A parameter fixed via set stays fixed in chained fits."""
    group = GroupFit(make_sweep())
    model = [(ModelName.Gaussian, dict(guess=True, set={"sigma": dict(value=0.2, vary=False)}))]
    results = group.fit(model, chain=True)
    assert all(r.fit_result.peak.values["sigma"] == 0.2 for r in results)


def test_table_flattens_columns_and_parameters() -> None:
    """Table gives one row per scan for each column, parameter and their errors."""
    group = GroupFit(make_sweep())
    group.fit(PEAK_AND_BACKGROUND)
    table = group.table(["temp"], ["center", "fwhm"])
    assert set(table) == {"uuid", "temp", "temp_err", "center", "center_err", "fwhm", "fwhm_err"}
    assert list(table["uuid"]) == ["s-10", "s-20", "s-30"]
    assert all(len(v) == 3 for v in table.values())


def test_summarize_before_fit_raises() -> None:
    """Summarizing with no results is an error, not an empty summary."""
    with pytest.raises(ValueError, match="fit"):
        GroupFit(make_sweep()).summarize("temp", "center")


def test_unknown_column_and_parameter_raise() -> None:
    """Misspelled names raise with the available choices."""
    group = GroupFit(make_sweep())
    group.fit([(ModelName.Gaussian, dict(guess=True))])
    with pytest.raises(KeyError, match="available"):
        group.summarize("field", "center")
    with pytest.raises(KeyError, match="available"):
        group.summarize("temp", "gamma")
