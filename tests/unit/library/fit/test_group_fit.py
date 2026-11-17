"""Tests for GroupFit: per-column fits and averages over a group of scans."""

import numpy as np
import pytest

from tavi.library.data.scan import UUID, Provenance, RawScan, ScanData, ScanMetadata, TaviMetadata
from tavi.library.fit import ColumnAverage, FitResult, GroupFit, GroupFitParams, GroupFitResult, ModelName

S1 = np.linspace(-2.0, 3.0, 101)


def make_scan(
    uuid_val: str,
    temp: float = 10.0,
    center: float = 0.0,
    sigma: float = 0.2,
    frequency: float = 2.0,
    seed: int = 0,
) -> RawScan:
    """A scan with a Poisson Gaussian ``detector``, a Poisson cosine ``pol`` and a near-constant ``temp``, vs s1."""
    rng = np.random.default_rng(seed)
    peak = 10.0 + 1000.0 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-((S1 - center) ** 2) / (2 * sigma**2))
    # Real counts are noisy; noiseless data makes lmfit's covariance degenerate (no stderr).
    detector = rng.poisson(peak).astype(float)
    pol = rng.poisson(500.0 + 300.0 * np.cos(frequency * (S1 + 0.3))).astype(float)
    return RawScan(
        uuid=UUID(value=uuid_val),
        data=ScanData(
            data={
                "s1": S1.tolist(),
                "s1_doubled": (2 * S1).tolist(),
                "detector": detector.tolist(),
                "pol": pol.tolist(),
                "temp": (temp + 0.01 * np.sin(S1)).tolist(),
                "_2theta": [40.0] * len(S1),
            }
        ),
        metadata=ScanMetadata(),
        tavimeta=TaviMetadata(default_axis=("s1", "detector"), friendly_name=uuid_val, friendly_path="/exp1"),
        prov=Provenance(raw_file=f"{uuid_val}.dat", contributing_scans={UUID(value=uuid_val): 1}),
    )


def make_sweep() -> list[RawScan]:
    """Three scans whose peak walks right and whose cosine speeds up as temperature rises."""
    return [
        make_scan("s-10", temp=10.0, center=0.0, frequency=2.0, seed=0),
        make_scan("s-20", temp=20.0, center=0.3, frequency=2.5, seed=1),
        make_scan("s-30", temp=30.0, center=0.6, frequency=3.0, seed=2),
    ]


def make_drift() -> tuple[list[RawScan], list[float]]:
    """Eight scans of a narrow peak drifting ~12 widths away from where an unchained fit starts."""
    centers = [2.0 + 0.05 * i for i in range(8)]
    return [make_scan(f"d{i}", center=c, sigma=0.03, seed=i) for i, c in enumerate(centers)], centers


PEAK_AND_BACKGROUND = [
    (ModelName.Gaussian, dict(guess=True, prefix="p_")),
    (ModelName.Linear, dict(slope=0.0, intercept=10.0, prefix="bg_")),
]
COSINE = [(ModelName.Cosine, dict(guess=True))]
DRIFT_START = [(ModelName.Gaussian, dict(center=2.0, sigma=0.03, amplitude=1000.0))]


def test_fit_returns_one_result_per_scan_in_order() -> None:
    """Each scan gets a GroupFitResult with its uuid, a FitResult per fitted column and a ColumnAverage per average."""
    scans = make_sweep()
    params = GroupFitParams(fits={"detector": PEAK_AND_BACKGROUND, "pol": COSINE}, averages=["temp"])
    results = GroupFit(scans, params).fit()
    assert [r.uuid for r in results] == [s.uuid for s in scans]
    for result in results:
        assert isinstance(result, GroupFitResult)
        assert set(result.fits) == {"detector", "pol"}
        assert all(isinstance(fit, FitResult) for fit in result.fits.values())
        assert isinstance(result.averages["temp"], ColumnAverage)


def test_peak_and_cosine_columns_fit_against_one_x() -> None:
    """A peak column and a cosine column of the same scans each recover their injected parameters."""
    params = GroupFitParams(fits={"detector": PEAK_AND_BACKGROUND, "pol": COSINE})
    results = GroupFit(make_sweep(), params).fit()
    assert np.allclose([r.fits["detector"]["p_"]["center"] for r in results], [0.0, 0.3, 0.6], atol=0.01)
    assert np.allclose([r.fits["pol"].frequency for r in results], [2.0, 2.5, 3.0], atol=0.05)
    assert np.allclose([r.fits["pol"].amplitude for r in results], 300.0, rtol=0.1)


def test_averages_give_mean_and_spread() -> None:
    """An averaged column is reduced to its per-scan mean and standard deviation."""
    results = GroupFit(make_sweep(), GroupFitParams(averages=["temp"])).fit()
    assert np.allclose([r.averages["temp"].mean for r in results], [10.0, 20.0, 30.0], atol=0.01)
    assert all(r.averages["temp"].std > 0 for r in results)
    assert all(r.fits == {} for r in results)


def test_x_col_defaults_to_scan_default_axis() -> None:
    """With no x_col the scan's default x axis is used; an explicit x_col is honored."""
    scans = make_sweep()
    default = GroupFit(scans, GroupFitParams(fits={"detector": PEAK_AND_BACKGROUND})).fit()
    doubled = GroupFit(scans, GroupFitParams(x_col="s1_doubled", fits={"detector": PEAK_AND_BACKGROUND})).fit()
    centers = np.array([r.fits["detector"]["p_"]["center"] for r in default])
    assert np.allclose([r.fits["detector"]["p_"]["center"] for r in doubled], 2 * centers, atol=0.02)


def test_digit_leading_column_accepts_underscore_spelling() -> None:
    """A column stored as ``_2theta`` can be asked for as ``2theta``."""
    results = GroupFit(make_sweep(), GroupFitParams(averages=["2theta"])).fit()
    assert results[0].averages["2theta"].mean == 40.0


def test_empty_params_raise() -> None:
    """Nothing to fit or average is an error, not an empty result."""
    with pytest.raises(ValueError, match="no columns"):
        GroupFit(make_sweep(), GroupFitParams()).fit()


def test_unknown_column_raises_with_available() -> None:
    """A misspelled fitted, averaged or x column raises listing the scan's columns."""
    for params in (
        GroupFitParams(fits={"detecter": COSINE}),
        GroupFitParams(averages=["field"]),
        GroupFitParams(x_col="s2", fits={"pol": COSINE}),
    ):
        with pytest.raises(KeyError, match="available"):
            GroupFit(make_sweep(), params).fit()


def test_chain_tracks_a_drifting_peak() -> None:
    """Chaining seeds each scan from its neighbour's fit, so a narrow drifting peak is tracked; unchained it is lost."""
    scans, centers = make_drift()
    unchained = GroupFit(scans, GroupFitParams(fits={"detector": DRIFT_START}, chain=False)).fit()
    assert not np.allclose([r.fits["detector"].center for r in unchained], centers, atol=0.01)
    chained = GroupFit(scans, GroupFitParams(fits={"detector": DRIFT_START})).fit()
    assert np.allclose([r.fits["detector"].center for r in chained], centers, atol=0.01)


def test_chain_keeps_set_options() -> None:
    """A parameter fixed via set stays fixed in every chained fit."""
    model = [(ModelName.Gaussian, dict(guess=True, set={"sigma": dict(value=0.2, vary=False)}))]
    results = GroupFit(make_sweep(), GroupFitParams(fits={"detector": model})).fit()
    assert all(r.fits["detector"].sigma == 0.2 for r in results)


def test_model_for_shows_chained_seed() -> None:
    """The first scan uses the group model; a later scan starts from its neighbour's fitted values, guess dropped."""
    scans = make_sweep()
    group = GroupFit(scans, GroupFitParams(fits={"detector": PEAK_AND_BACKGROUND}))
    results = group.fit()
    assert group.model_for(scans[0].uuid, "detector") == PEAK_AND_BACKGROUND
    (_, peak), (_, background) = group.model_for(scans[1].uuid, "detector")
    assert "guess" not in peak
    assert peak["center"] == results[0].fits["detector"]["p_"]["center"]
    assert set(peak) == {"prefix", "center", "sigma", "amplitude"}
    assert background["slope"] == results[0].fits["detector"]["bg_"]["slope"]


def test_non_finite_seed_falls_back_to_group_model() -> None:
    """A NaN in the neighbour's fit is not seeded; the component keeps its own start."""
    scans = make_sweep()
    group = GroupFit(scans, GroupFitParams(fits={"detector": PEAK_AND_BACKGROUND}))
    group.fit()
    group.results[0].fits["detector"]["p_"].values["center"] = np.nan
    (_, peak), _ = group.model_for(scans[1].uuid, "detector")
    assert peak == PEAK_AND_BACKGROUND[0][1]


def test_override_applies_to_one_scan_only() -> None:
    """An override fixes a parameter in its own scan; the other scans fit freely."""
    scans = make_sweep()
    group = GroupFit(scans, GroupFitParams(fits={"detector": PEAK_AND_BACKGROUND}, chain=False))
    group.override(scans[1].uuid, "detector", "p_", {"set": {"sigma": {"value": 0.15, "vary": False}}})
    results = group.fit()
    sigmas = [r.fits["detector"]["p_"]["sigma"] for r in results]
    assert sigmas[1] == 0.15
    assert sigmas[0] != 0.15 and sigmas[2] != 0.15


def test_up_front_overrides_match_override_calls() -> None:
    """Overrides given in GroupFitParams behave like override() calls, without mutating the caller's params."""
    scans = make_sweep()
    patch = {"set": {"sigma": {"value": 0.15, "vary": False}}}
    params = GroupFitParams(
        fits={"detector": PEAK_AND_BACKGROUND}, overrides={scans[2].uuid: {"detector": {"p_": patch}}}
    )
    up_front = GroupFit(scans, params)
    later = GroupFit(scans, GroupFitParams(fits={"detector": PEAK_AND_BACKGROUND}))
    later.override(scans[2].uuid, "detector", "p_", patch)
    assert up_front.model_for(scans[2].uuid, "detector") == later.model_for(scans[2].uuid, "detector")

    up_front.override(scans[0].uuid, "detector", "p_", patch)
    assert set(params.overrides) == {scans[2].uuid}


def test_overrides_merge_set_per_parameter() -> None:
    """Successive overrides accumulate, and a set option on one parameter keeps that parameter's other options."""
    scans = make_sweep()
    model = [(ModelName.Gaussian, dict(guess=True, set={"sigma": {"min": 0.0}}))]
    group = GroupFit(scans, GroupFitParams(fits={"detector": model}, chain=False))
    group.override(scans[0].uuid, "detector", "", {"set": {"sigma": {"vary": False}}})
    group.override(scans[0].uuid, "detector", "", {"set": {"center": {"max": 0.5}}})
    ((_, params),) = group.model_for(scans[0].uuid, "detector")
    assert params["set"] == {"sigma": {"min": 0.0, "vary": False}, "center": {"max": 0.5}}
    assert model[0][1]["set"] == {"sigma": {"min": 0.0}}


def test_override_wins_over_chained_seed() -> None:
    """An override's starting value replaces the value seeded from the previous scan."""
    scans = make_sweep()
    group = GroupFit(scans, GroupFitParams(fits={"detector": PEAK_AND_BACKGROUND}))
    group.override(scans[1].uuid, "detector", "p_", {"center": 0.25})
    group.fit()
    (_, peak), _ = group.model_for(scans[1].uuid, "detector")
    assert peak["center"] == 0.25


def test_override_rejects_unknown_scan_column_or_prefix() -> None:
    """An override naming a scan, column or component the group does not have raises."""
    scans = make_sweep()
    group = GroupFit(scans, GroupFitParams(fits={"detector": PEAK_AND_BACKGROUND}))
    with pytest.raises(KeyError, match="not in this group"):
        group.override(UUID(value="nope"), "detector", "p_", {})
    with pytest.raises(KeyError):
        group.override(scans[0].uuid, "pol", "p_", {})
    with pytest.raises(KeyError, match="available"):
        group.override(scans[0].uuid, "detector", "q_", {})
    assert group.group_fit_params.overrides == {}


def test_up_front_override_is_validated() -> None:
    """A typo in GroupFitParams.overrides fails when the GroupFit is built, not silently at fit time."""
    params = GroupFitParams(
        fits={"detector": PEAK_AND_BACKGROUND}, overrides={UUID(value="s-10"): {"detector": {"x_": {}}}}
    )
    with pytest.raises(KeyError, match="available"):
        GroupFit(make_sweep(), params)
