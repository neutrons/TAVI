"""Tests for FitModel."""

import math

import numpy as np
import numpy.testing as npt
import pytest

from tavi.backend.model.fit_model import FitModel
from tavi.library.data.fit_entry import (
    FitEntry,
    FitMember,
    FitRequest,
    FitSpec,
    ParamField,
    PeakField,
    SuggestBackgroundParamsRequest,
    SuggestPeakParamsRequest,
)
from tavi.library.data.model_response import ResponseCode
from tavi.library.data.plot import PlotSeries
from tavi.library.data.scan import UUID, Provenance, RawScan, ScanData, ScanMetadata, TaviMetadata
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.exception_event import ExceptionEvent
from tavi.meta.event.type.presenter_event import (
    BackgroundParamsSuggestedEvent,
    FitRecomputeEvent,
    PeakParamsSuggestedEvent,
    SaveFitEvent,
    SyncFitEvent,
    SyncFitSpecEvent,
)

GAUSSIAN_AMPLITUDE = 5.0
GAUSSIAN_CENTER = 0.0
GAUSSIAN_SIGMA = 1.0
GAUSSIAN_FWHM = GAUSSIAN_SIGMA * 2 * math.sqrt(2 * math.log(2))
# lmfit's own "amplitude" is the integrated area (height * sigma * sqrt(2*pi)), not the peak
# height make_gaussian_xy is parametrized by - what FitResultSummary.amplitude reports.
GAUSSIAN_AREA = GAUSSIAN_AMPLITUDE * GAUSSIAN_SIGMA * math.sqrt(2 * math.pi)


def make_param(value, fixed=False, minimum="", maximum="") -> ParamField:
    return ParamField(value=str(value), fixed=fixed, minimum=str(minimum), maximum=str(maximum))


def make_gaussian_xy(n=101, center=GAUSSIAN_CENTER, amplitude=GAUSSIAN_AMPLITUDE):
    x = np.linspace(-5, 5, n)
    y = amplitude * np.exp(-((x - center) ** 2) / (2 * GAUSSIAN_SIGMA**2))
    return x, y


def make_scan(uuid_val="scan-001", x=None, y=None, friendly_name="test_scan") -> RawScan:
    """A scan whose ``qh``/``detector`` columns default to the noiseless reference Gaussian."""
    default_x, default_y = make_gaussian_xy()
    x = default_x if x is None else np.asarray(x)
    y = default_y if y is None else np.asarray(y)
    return RawScan(
        uuid=UUID(value=uuid_val),
        data=ScanData(data={"qh": x.tolist(), "detector": y.tolist()}),
        metadata=ScanMetadata(),
        tavimeta=TaviMetadata(default_axis=("qh", "detector"), friendly_name=friendly_name, friendly_path="/exp1"),
        prov=Provenance(raw_file="scan.dat", contributing_scans={UUID(value=uuid_val): 1}),
    )


def put_scans(raw_scans: dict, *scans: RawScan) -> None:
    """Load ``scans`` into the FitModel's live raw_scans handle, replacing any with the same uuid."""
    for scan in scans:
        raw_scans[scan.uuid] = scan


def make_series(uuid_val="scan-001", friendly_name=None) -> PlotSeries:
    return PlotSeries(
        source_scan_uuid=UUID(value=uuid_val),
        scan_name="test_scan",
        friendly_name=friendly_name,
        normalized_by=None,
        x_name="qh",
        y_name="detector",
        error_name="err",
    )


def make_spec(**overrides) -> FitSpec:
    defaults = dict(
        range_min="-5",
        range_max="5",
        background="None",
        background_constant=make_param(0),
        peaks=[
            PeakField(
                shape="Gaussian",
                amplitude=make_param(GAUSSIAN_AMPLITUDE),
                center=make_param(GAUSSIAN_CENTER),
                fwhm=make_param(GAUSSIAN_FWHM),
            )
        ],
    )
    defaults.update(overrides)
    return FitSpec(**defaults)


def make_request(series=None, seed_from_previous=False, fit_uuid=None, **spec_overrides) -> FitRequest:
    return FitRequest(
        spec=make_spec(**spec_overrides),
        series=[make_series()] if series is None else series,
        seed_from_previous=seed_from_previous,
        fit_uuid=fit_uuid,
    )


def make_member(series=None, **spec_overrides) -> FitMember:
    return FitMember(series=series or make_series(), **make_spec(**spec_overrides).model_dump())


def make_fit_entry(members=None, uuid_val="fit-001") -> FitEntry:
    return FitEntry(uuid=UUID(value=uuid_val), members=members or [make_member()])


@pytest.fixture
def raw_scans():
    scans = {}
    put_scans(scans, make_scan())
    return scans


@pytest.fixture
def fits():
    return {}


@pytest.fixture
def model(raw_scans, fits):
    return FitModel(raw_scans=raw_scans, fits=fits)


def test_perform_fit_returns_ok(model):
    response = model.perform_fit(make_request())
    assert response.code == ResponseCode.OK


def test_perform_fit_publishes_sync_fit_event(model):
    received = []
    EventBroker().register(SyncFitEvent, received.append)

    model.perform_fit(make_request())

    assert len(received) == 1
    assert len(received[0].outcomes) == 1


def test_perform_fit_does_not_cache_the_curve_on_the_member(model):
    """The fit's member must not carry the fitted (x, best_fit) points - only its own series/spec/result."""
    received = []
    EventBroker().register(SyncFitEvent, received.append)

    model.perform_fit(make_request())

    member = received[0].outcomes[0].member
    assert not hasattr(member, "x")
    assert not hasattr(member, "best_fit")
    assert member.source_scan_uuid == UUID(value="scan-001")


def test_perform_fit_curve_carries_its_source_scan(model):
    """The curve identifies the series it belongs to, so redrawing it replaces that series' line."""
    received = []
    EventBroker().register(SyncFitEvent, received.append)

    model.perform_fit(make_request())

    assert received[0].outcomes[0].curve.source_scan_uuid == UUID(value="scan-001")


def test_perform_fit_publishes_each_component_separately(model, raw_scans):
    """Every component of a composite model has to reach the curve - "Plot Separately" draws them."""
    received = []
    EventBroker().register(SyncFitEvent, received.append)

    model.perform_fit(make_linear_background_request(raw_scans))

    assert sorted(received[0].outcomes[0].curve.components) == ["bg_", "peak1_"]


def test_perform_fit_components_sum_to_the_composite_curve(model, raw_scans):
    received = []
    EventBroker().register(SyncFitEvent, received.append)

    model.perform_fit(make_linear_background_request(raw_scans))

    curve = received[0].outcomes[0].curve
    total = np.sum([np.array(values) for values in curve.components.values()], axis=0)
    npt.assert_allclose(total, curve.best_fit)


def test_perform_fit_components_are_evaluated_on_the_curves_own_grid(model, raw_scans):
    received = []
    EventBroker().register(SyncFitEvent, received.append)

    model.perform_fit(make_linear_background_request(raw_scans))

    curve = received[0].outcomes[0].curve
    assert all(len(values) == len(curve.x) for values in curve.components.values())


def test_perform_fit_single_component_publishes_no_components(model):
    """One component *is* the composite curve - separating it out would just overplot it."""
    received = []
    EventBroker().register(SyncFitEvent, received.append)

    model.perform_fit(make_request(background="None"))

    assert received[0].outcomes[0].curve.components == {}


def test_perform_fit_multi_peak_publishes_one_component_per_peak(model):
    received = []
    EventBroker().register(SyncFitEvent, received.append)

    peak = PeakField(shape="Gaussian", amplitude=make_param(""), center=make_param(0), fwhm=make_param(1))
    model.perform_fit(make_request(background="None", peaks=[peak, peak.model_copy(deep=True)]))

    assert sorted(received[0].outcomes[0].curve.components) == ["peak1_", "peak2_"]


def test_perform_fit_publishes_the_data_it_was_fit_against(model):
    received = []
    EventBroker().register(SyncFitEvent, received.append)

    model.perform_fit(make_request())

    x, y = make_gaussian_xy()
    data = received[0].outcomes[0].data
    npt.assert_allclose(data.x, x)
    npt.assert_allclose(data.y, y)
    npt.assert_allclose(data.err, np.sqrt(np.abs(y)))


def test_perform_fit_without_a_fit_uuid_mints_a_new_one_each_time(model):
    received = []
    EventBroker().register(SyncFitEvent, received.append)

    model.perform_fit(make_request())
    model.perform_fit(make_request())

    assert received[0].fit_uuid != received[1].fit_uuid


def test_perform_fit_with_a_fit_uuid_refits_that_fit_in_place(model):
    """A refit carries the existing uuid through, so TaviProjectModel overwrites that entry."""
    received = []
    EventBroker().register(SyncFitEvent, received.append)

    model.perform_fit(make_request(fit_uuid=UUID(value="fit-007")))

    assert received[0].fit_uuid == UUID(value="fit-007")


def test_perform_fit_recovers_amplitude(model):
    received = []
    EventBroker().register(SyncFitEvent, received.append)

    model.perform_fit(make_request())

    curve = received[0].outcomes[0].curve
    assert max(curve.best_fit) == pytest.approx(GAUSSIAN_AMPLITUDE, rel=1e-3)


def test_perform_fit_reduced_chi_squared_near_zero_for_noiseless_data(model):
    received = []
    EventBroker().register(SyncFitEvent, received.append)

    model.perform_fit(make_request())

    assert received[0].outcomes[0].member.result.reduced_chi_squared == pytest.approx(0.0, abs=1e-6)


def test_perform_fit_result_recovers_center_and_fwhm(model):
    received = []
    EventBroker().register(SyncFitEvent, received.append)

    model.perform_fit(make_request())

    result = received[0].outcomes[0].member.result
    assert result.peaks[0].amplitude == pytest.approx(GAUSSIAN_AREA, rel=1e-3)
    assert result.peaks[0].center == pytest.approx(GAUSSIAN_CENTER, abs=1e-3)
    assert result.peaks[0].fwhm == pytest.approx(GAUSSIAN_FWHM, rel=1e-3)


def test_perform_fit_reports_uncertainties_for_noisy_data(model, raw_scans):
    """
    Perfectly noiseless synthetic data can leave lmfit's covariance degenerate (stderr=None) -
    add a touch of deterministic noise so there's something real to estimate an uncertainty from.
    """
    received = []
    EventBroker().register(SyncFitEvent, received.append)

    _x, y = make_gaussian_xy()
    rng = np.random.default_rng(seed=0)
    put_scans(raw_scans, make_scan(y=y + rng.normal(scale=0.02, size=y.shape)))
    model.perform_fit(make_request())

    result = received[0].outcomes[0].member.result
    assert result.peaks[0].amplitude_err is not None
    assert result.peaks[0].center_err is not None
    assert result.peaks[0].fwhm_err is not None


def test_perform_fit_flat_background_is_linear_with_its_slope_fixed(model, raw_scans):
    """There is no separate "Constant" model - a flat background is Linear with its slope fixed at 0."""
    received = []
    EventBroker().register(SyncFitEvent, received.append)

    _x, y = make_gaussian_xy()
    put_scans(raw_scans, make_scan(y=y + 2.0))
    request = make_request(
        background="Linear",
        background_constant=make_param(0),
        background_slope=make_param(0, fixed=True),
    )

    model.perform_fit(request)

    assert len(received) == 1
    result = received[0].outcomes[0].member.result
    assert result.background_constant == pytest.approx(2.0, abs=1e-2)
    assert result.background_slope == pytest.approx(0.0, abs=1e-9)


def test_perform_fit_trims_to_range(model):
    received = []
    EventBroker().register(SyncFitEvent, received.append)

    model.perform_fit(make_request(range_min="-1", range_max="1"))

    curve = received[0].outcomes[0].curve
    assert min(curve.x) >= -1
    assert max(curve.x) <= 1
    # The curve is re-evaluated on a fine grid rather than at the scan's own points, so its
    # point count says nothing about trimming - its span does.
    assert max(curve.x) - min(curve.x) == pytest.approx(2.0, abs=0.1)


def test_perform_fit_unsupported_peak_shape_reports_error(model):
    errors = []
    EventBroker().register(ExceptionEvent, errors.append)
    computed = []
    EventBroker().register(SyncFitEvent, computed.append)

    request = make_request(
        peaks=[PeakField(shape="Pseudo-Voigt", amplitude=make_param(1), center=make_param(0), fwhm=make_param(1))]
    )
    model.perform_fit(request)

    assert len(errors) == 1
    assert computed == []


def test_perform_fit_unsupported_background_reports_error(model):
    errors = []
    EventBroker().register(ExceptionEvent, errors.append)
    computed = []
    EventBroker().register(SyncFitEvent, computed.append)

    model.perform_fit(make_request(background="Quadratic"))

    assert len(errors) == 1
    assert computed == []


def test_perform_fit_non_numeric_peak_field_reports_error(model):
    errors = []
    EventBroker().register(ExceptionEvent, errors.append)

    request = make_request(
        peaks=[
            PeakField(shape="Gaussian", amplitude=make_param("not-a-number"), center=make_param(0), fwhm=make_param(1))
        ]
    )
    model.perform_fit(request)

    assert len(errors) == 1


def test_perform_fit_non_numeric_range_reports_error(model):
    errors = []
    EventBroker().register(ExceptionEvent, errors.append)

    model.perform_fit(make_request(range_min="abc"))

    assert len(errors) == 1


def test_perform_fit_range_with_too_few_points_reports_error(model):
    errors = []
    EventBroker().register(ExceptionEvent, errors.append)

    model.perform_fit(make_request(range_min="4.99", range_max="5.0"))

    assert len(errors) == 1


def test_perform_fit_unloaded_scan_reports_error(model):
    errors = []
    EventBroker().register(ExceptionEvent, errors.append)
    computed = []
    EventBroker().register(SyncFitEvent, computed.append)

    model.perform_fit(make_request(series=[make_series("scan-missing")]))

    assert len(errors) == 1
    assert computed == []


# ---------------------------------------------------------------------------
# Blank initial params -> auto-guess from data
# ---------------------------------------------------------------------------


def test_perform_fit_guesses_blank_peak_params(model):
    """Leaving all three peak fields blank must guess from the data, not error out."""
    received = []
    EventBroker().register(SyncFitEvent, received.append)
    errors = []
    EventBroker().register(ExceptionEvent, errors.append)

    request = make_request(
        peaks=[PeakField(shape="Gaussian", amplitude=make_param(""), center=make_param(""), fwhm=make_param(""))]
    )
    model.perform_fit(request)

    assert errors == []
    assert len(received) == 1
    result = received[0].outcomes[0].member.result
    assert result.peaks[0].amplitude == pytest.approx(GAUSSIAN_AREA, rel=1e-2)
    assert result.peaks[0].center == pytest.approx(GAUSSIAN_CENTER, abs=1e-2)


def test_perform_fit_guesses_only_the_blank_peak_field(model):
    """A deliberately-wrong explicit center must survive - only the blank fields get guessed."""
    received = []
    EventBroker().register(SyncFitEvent, received.append)

    request = make_request(
        peaks=[
            PeakField(
                shape="Gaussian", amplitude=make_param(""), center=make_param(GAUSSIAN_CENTER), fwhm=make_param("")
            )
        ]
    )
    model.perform_fit(request)

    assert received[0].outcomes[0].member.result.peaks[0].center == pytest.approx(GAUSSIAN_CENTER, abs=1e-3)


def test_perform_fit_guesses_blank_background_constant(model, raw_scans):
    received = []
    EventBroker().register(SyncFitEvent, received.append)

    _x, y = make_gaussian_xy()
    put_scans(raw_scans, make_scan(y=y + 2.0))
    request = make_request(
        background="Linear",
        background_constant=make_param(""),
        background_slope=make_param(0, fixed=True),
    )
    model.perform_fit(request)

    assert received[0].outcomes[0].member.result.background_constant == pytest.approx(2.0, abs=0.5)


# ---------------------------------------------------------------------------
# Sequential fit - several series in one request
# ---------------------------------------------------------------------------

SEQUENTIAL_CENTERS = [0.25, 0.75, 1.25]


def make_sequential_scans(raw_scans, count=3, centers=SEQUENTIAL_CENTERS):
    """Load ``count`` scans whose peak drifts a little per scan, and return a series for each, in order."""
    series = []
    for index in range(count):
        uuid_val = f"scan-{index:03d}"
        x, y = make_gaussian_xy(center=centers[index])
        put_scans(raw_scans, make_scan(uuid_val, x=x, y=y, friendly_name=f"run-{index}"))
        series.append(make_series(uuid_val, friendly_name=f"run-{index}"))
    return series


def test_sequential_fit_publishes_one_event_with_one_outcome_per_series_in_order(model, raw_scans):
    received = []
    EventBroker().register(SyncFitEvent, received.append)
    series = make_sequential_scans(raw_scans)

    model.perform_fit(make_request(series=series, seed_from_previous=True))

    assert len(received) == 1
    assert [outcome.member.source_scan_uuid for outcome in received[0].outcomes] == [s.source_scan_uuid for s in series]
    assert all(outcome.member.result is not None for outcome in received[0].outcomes)


def test_sequential_fit_fits_each_series_against_its_own_data(model, raw_scans):
    received = []
    EventBroker().register(SyncFitEvent, received.append)
    series = make_sequential_scans(raw_scans)

    model.perform_fit(make_request(series=series, seed_from_previous=True))

    centers = [outcome.member.result.peaks[0].center for outcome in received[0].outcomes]
    assert centers == pytest.approx(SEQUENTIAL_CENTERS, abs=1e-3)


def test_sequential_fit_does_not_freeze_a_parameter_that_fits_to_zero(model, raw_scans):
    # The first center fits to ~1e-13; seeding the next member with that verbatim gives MINPACK a
    # vanishing relative step and leaves every later center stuck there.
    centers = [0.0, 0.5, 1.0]
    received = []
    EventBroker().register(SyncFitEvent, received.append)
    series = make_sequential_scans(raw_scans, centers=centers)

    model.perform_fit(make_request(series=series, seed_from_previous=True))

    fitted = [outcome.member.result.peaks[0].center for outcome in received[0].outcomes]
    assert fitted == pytest.approx(centers, abs=1e-3)
    assert received[0].outcomes[1].member.peaks[0].center.value == "0.0"


def test_sequential_fit_seeds_each_member_from_the_previous_result(model, raw_scans):
    received = []
    EventBroker().register(SyncFitEvent, received.append)
    series = make_sequential_scans(raw_scans, count=2)

    model.perform_fit(make_request(series=series, seed_from_previous=True))

    first, second = received[0].outcomes
    assert second.member.peaks[0].center.value == repr(first.member.result.peaks[0].center)
    assert second.member.peaks[0].amplitude.value == repr(first.member.result.peaks[0].amplitude)
    assert second.member.peaks[0].fwhm.value == repr(first.member.result.peaks[0].fwhm)


def test_sequential_fit_first_member_keeps_the_requested_spec(model, raw_scans):
    received = []
    EventBroker().register(SyncFitEvent, received.append)
    series = make_sequential_scans(raw_scans, count=2)

    model.perform_fit(make_request(series=series, seed_from_previous=True))

    assert received[0].outcomes[0].member.peaks[0].center.value == str(GAUSSIAN_CENTER)


def test_fit_without_seeding_applies_the_spec_as_is_to_every_series(model, raw_scans):
    received = []
    EventBroker().register(SyncFitEvent, received.append)
    series = make_sequential_scans(raw_scans, count=2)

    model.perform_fit(make_request(series=series, seed_from_previous=False))

    assert [outcome.member.peaks[0].center.value for outcome in received[0].outcomes] == [str(GAUSSIAN_CENTER)] * 2


def test_sequential_fit_seeding_keeps_fixed_flags_and_bounds(model, raw_scans):
    """Only the values carry over - fixed flags and bounds stay as the spec has them."""
    received = []
    EventBroker().register(SyncFitEvent, received.append)
    series = make_sequential_scans(raw_scans, count=2)
    peak = PeakField(
        shape="Gaussian",
        amplitude=make_param(GAUSSIAN_AMPLITUDE, minimum=0),
        center=make_param(GAUSSIAN_CENTER, minimum=-2, maximum=2),
        fwhm=make_param(GAUSSIAN_FWHM, fixed=True),
    )

    model.perform_fit(make_request(series=series, seed_from_previous=True, peaks=[peak]))

    seeded = received[0].outcomes[1].member.peaks[0]
    assert seeded.fwhm.fixed is True
    assert seeded.center.fixed is False
    assert (seeded.center.minimum, seeded.center.maximum) == ("-2", "2")
    assert (seeded.amplitude.minimum, seeded.amplitude.maximum) == ("0", "")
    assert received[0].outcomes[1].member.result.peaks[0].fwhm == pytest.approx(GAUSSIAN_FWHM, rel=1e-6)


def test_sequential_fit_failed_member_is_published_without_a_result(model, raw_scans):
    received = []
    EventBroker().register(SyncFitEvent, received.append)
    series = make_sequential_scans(raw_scans)
    series[1] = make_series("scan-missing", friendly_name="run-missing")

    model.perform_fit(make_request(series=series, seed_from_previous=True))

    outcomes = received[0].outcomes
    assert len(outcomes) == 3
    assert outcomes[1].member.source_scan_uuid == UUID(value="scan-missing")
    assert outcomes[1].member.result is None
    assert outcomes[1].curve is None
    assert outcomes[1].data is None


def test_sequential_fit_seeding_continues_from_the_last_good_member(model, raw_scans):
    received = []
    EventBroker().register(SyncFitEvent, received.append)
    series = make_sequential_scans(raw_scans)
    series[1] = make_series("scan-missing", friendly_name="run-missing")

    model.perform_fit(make_request(series=series, seed_from_previous=True))

    first, _failed, third = received[0].outcomes
    assert third.member.result is not None
    assert third.member.peaks[0].center.value == repr(first.member.result.peaks[0].center)


def test_sequential_fit_with_every_member_failing_publishes_nothing(model):
    computed = []
    EventBroker().register(SyncFitEvent, computed.append)
    series = [make_series("scan-missing-a", "run-a"), make_series("scan-missing-b", "run-b")]

    model.perform_fit(make_request(series=series, seed_from_previous=True))

    assert computed == []


def test_multi_series_errors_aggregate_into_one_report_prefixed_by_run_name(model, raw_scans):
    errors = []
    EventBroker().register(ExceptionEvent, errors.append)
    series = make_sequential_scans(raw_scans, count=2)

    model.perform_fit(make_request(series=series, range_min="abc"))

    assert len(errors) == 1
    lines = errors[0].error.message.splitlines()
    assert len(lines) == 2
    assert lines[0].startswith("run-0: ")
    assert lines[1].startswith("run-1: ")


def test_multi_series_fit_reuses_a_supplied_fit_uuid(model, raw_scans):
    received = []
    EventBroker().register(SyncFitEvent, received.append)
    series = make_sequential_scans(raw_scans, count=2)

    model.perform_fit(make_request(series=series, seed_from_previous=True, fit_uuid=UUID(value="fit-007")))

    assert received[0].fit_uuid == UUID(value="fit-007")


# ---------------------------------------------------------------------------
# FitRecomputeEvent - recompute a selected fit against live data, never a cached curve
# ---------------------------------------------------------------------------


def test_fit_recompute_recomputes_against_current_raw_scan_data(model):
    entry = make_fit_entry()

    received = []
    EventBroker().register(SyncFitEvent, received.append)
    EventBroker().publish(FitRecomputeEvent(fits=[entry]))

    assert len(received) == 1
    assert received[0].fit_uuid == entry.uuid
    assert received[0].outcomes[0].member.result.peaks[0].amplitude == pytest.approx(GAUSSIAN_AREA, rel=1e-2)


def test_fit_recompute_reflects_changed_underlying_data(model, raw_scans):
    """Recompute must reflect the scan's *current* data, not whatever it was fit against originally."""
    entry = make_fit_entry()

    _x, doubled = make_gaussian_xy(amplitude=2 * GAUSSIAN_AMPLITUDE)
    raw_scans[UUID(value="scan-001")].data.data["detector"] = doubled.tolist()

    received = []
    EventBroker().register(SyncFitEvent, received.append)
    EventBroker().publish(FitRecomputeEvent(fits=[entry]))

    assert received[0].outcomes[0].member.result.peaks[0].amplitude == pytest.approx(2 * GAUSSIAN_AREA, rel=1e-2)


def test_fit_recompute_unknown_source_scan_is_noop(model):
    entry = make_fit_entry(members=[make_member(make_series("scan-missing"))])

    received = []
    EventBroker().register(SyncFitEvent, received.append)
    EventBroker().publish(FitRecomputeEvent(fits=[entry]))

    assert received == []


def test_fit_recompute_publishes_every_member_of_a_multi_member_fit(model, raw_scans):
    series = make_sequential_scans(raw_scans)
    entry = make_fit_entry(members=[make_member(s) for s in series])

    received = []
    EventBroker().register(SyncFitEvent, received.append)
    EventBroker().publish(FitRecomputeEvent(fits=[entry]))

    assert len(received) == 1
    assert [outcome.member.source_scan_uuid for outcome in received[0].outcomes] == [s.source_scan_uuid for s in series]
    centers = [outcome.member.result.peaks[0].center for outcome in received[0].outcomes]
    assert centers == pytest.approx(SEQUENTIAL_CENTERS, abs=1e-3)


def test_fit_recompute_uses_each_members_own_spec_without_reseeding(model, raw_scans):
    series = make_sequential_scans(raw_scans, count=2)
    members = [make_member(series[0]), make_member(series[1], range_min="-1", range_max="1")]
    entry = make_fit_entry(members=members)

    received = []
    EventBroker().register(SyncFitEvent, received.append)
    EventBroker().publish(FitRecomputeEvent(fits=[entry]))

    first, second = received[0].outcomes
    assert first.member.peaks[0].center.value == str(GAUSSIAN_CENTER)
    assert second.member.peaks[0].center.value == str(GAUSSIAN_CENTER)
    assert max(second.curve.x) - min(second.curve.x) == pytest.approx(2.0, abs=0.1)


def test_fit_recompute_skips_members_whose_scan_is_not_loaded(model, raw_scans):
    series = make_sequential_scans(raw_scans, count=2)
    members = [make_member(series[0]), make_member(make_series("scan-missing")), make_member(series[1])]
    entry = make_fit_entry(members=members)

    received = []
    EventBroker().register(SyncFitEvent, received.append)
    EventBroker().publish(FitRecomputeEvent(fits=[entry]))

    assert [outcome.member.source_scan_uuid for outcome in received[0].outcomes] == [s.source_scan_uuid for s in series]


def test_fit_recompute_publishes_one_event_per_fit(model, raw_scans):
    series = make_sequential_scans(raw_scans, count=2)
    entries = [
        make_fit_entry(members=[make_member(series[0])], uuid_val="fit-a"),
        make_fit_entry(members=[make_member(series[1])], uuid_val="fit-b"),
    ]

    received = []
    EventBroker().register(SyncFitEvent, received.append)
    EventBroker().publish(FitRecomputeEvent(fits=entries))

    assert [event.fit_uuid for event in received] == [UUID(value="fit-a"), UUID(value="fit-b")]


# ---------------------------------------------------------------------------
# Linear background
# ---------------------------------------------------------------------------

LINEAR_BG_SLOPE = 0.4
LINEAR_BG_INTERCEPT = 2.0


def make_gaussian_on_a_slope_xy(n=201):
    """A Gaussian riding on a sloped baseline, so both components have a known true answer."""
    x = np.linspace(-5, 5, n)
    peak = GAUSSIAN_AMPLITUDE * np.exp(-((x - GAUSSIAN_CENTER) ** 2) / (2 * GAUSSIAN_SIGMA**2))
    return x, peak + LINEAR_BG_SLOPE * x + LINEAR_BG_INTERCEPT


def make_linear_background_request(raw_scans, **overrides) -> FitRequest:
    """Load the sloped-baseline scan in place of the reference one and return a Linear-background request."""
    x, y = make_gaussian_on_a_slope_xy()
    put_scans(raw_scans, make_scan(x=x, y=y))
    defaults = dict(
        background="Linear",
        background_constant=make_param(""),
        background_slope=make_param(""),
    )
    defaults.update(overrides)
    return make_request(**defaults)


def test_perform_fit_linear_background_is_supported(model, raw_scans):
    errors = []
    EventBroker().register(ExceptionEvent, errors.append)
    computed = []
    EventBroker().register(SyncFitEvent, computed.append)

    model.perform_fit(make_linear_background_request(raw_scans))

    assert errors == []
    assert len(computed) == 1


def test_perform_fit_linear_background_recovers_slope_and_intercept(model, raw_scans):
    computed = []
    EventBroker().register(SyncFitEvent, computed.append)

    model.perform_fit(make_linear_background_request(raw_scans))

    result = computed[0].outcomes[0].member.result
    assert result.background_slope == pytest.approx(LINEAR_BG_SLOPE, abs=0.05)
    assert result.background_constant == pytest.approx(LINEAR_BG_INTERCEPT, abs=0.1)


def test_perform_fit_linear_background_still_recovers_the_peak(model, raw_scans):
    computed = []
    EventBroker().register(SyncFitEvent, computed.append)

    model.perform_fit(make_linear_background_request(raw_scans))

    result = computed[0].outcomes[0].member.result
    assert result.peaks[0].center == pytest.approx(GAUSSIAN_CENTER, abs=0.05)
    assert result.peaks[0].amplitude == pytest.approx(GAUSSIAN_AREA, rel=0.05)


def test_perform_fit_without_a_background_reports_no_background_terms(model):
    """A "None" background has nothing to report - both terms stay unset rather than reading as a fitted 0."""
    computed = []
    EventBroker().register(SyncFitEvent, computed.append)

    model.perform_fit(make_request(background="None"))

    result = computed[0].outcomes[0].member.result
    assert result.background_constant is None
    assert result.background_slope is None


def test_perform_fit_linear_background_honours_a_fixed_slope(model, raw_scans):
    computed = []
    EventBroker().register(SyncFitEvent, computed.append)

    model.perform_fit(make_linear_background_request(raw_scans, background_slope=make_param(0, fixed=True)))

    assert computed[0].outcomes[0].member.result.background_slope == pytest.approx(0.0, abs=1e-9)


def test_perform_fit_linear_background_non_numeric_slope_reports_error(model, raw_scans):
    errors = []
    EventBroker().register(ExceptionEvent, errors.append)
    computed = []
    EventBroker().register(SyncFitEvent, computed.append)

    model.perform_fit(make_linear_background_request(raw_scans, background_slope=make_param("abc")))

    assert len(errors) == 1
    assert computed == []


# ---------------------------------------------------------------------------
# suggest_peak_params
# ---------------------------------------------------------------------------


def make_suggest_request(**overrides):
    x, y = make_gaussian_xy()
    defaults = dict(
        source_scan_uuid=UUID(value="scan-001"),
        x=x.tolist(),
        y=y.tolist(),
        range_min="-5",
        range_max="5",
        shape="Gaussian",
    )
    defaults.update(overrides)
    return SuggestPeakParamsRequest(**defaults)


def test_suggest_peak_params_returns_ok(model):
    response = model.suggest_peak_params(make_suggest_request())
    assert response.code == ResponseCode.OK


def test_suggest_peak_params_publishes_event(model):
    received = []
    EventBroker().register(PeakParamsSuggestedEvent, received.append)

    model.suggest_peak_params(make_suggest_request())

    assert len(received) == 1
    assert received[0].source_scan_uuid == UUID(value="scan-001")


def test_suggest_peak_params_finds_center_near_true_peak(model):
    received = []
    EventBroker().register(PeakParamsSuggestedEvent, received.append)

    model.suggest_peak_params(make_suggest_request())

    assert received[0].center == pytest.approx(GAUSSIAN_CENTER, abs=0.1)


def test_suggest_peak_params_amplitude_is_positive(model):
    received = []
    EventBroker().register(PeakParamsSuggestedEvent, received.append)

    model.suggest_peak_params(make_suggest_request())

    assert received[0].amplitude > 0


def test_suggest_peak_params_fwhm_is_positive(model):
    received = []
    EventBroker().register(PeakParamsSuggestedEvent, received.append)

    model.suggest_peak_params(make_suggest_request())

    assert received[0].fwhm > 0


def test_suggest_peak_params_works_for_lorentzian_and_voigt(model):
    for shape in ("Lorentzian", "Voigt"):
        received = []
        EventBroker().register(PeakParamsSuggestedEvent, received.append)

        model.suggest_peak_params(make_suggest_request(shape=shape))

        assert len(received) == 1
        EventBroker().registry[PeakParamsSuggestedEvent].remove(received.append)


def test_suggest_peak_params_unsupported_shape_reports_error(model):
    errors = []
    EventBroker().register(ExceptionEvent, errors.append)
    suggested = []
    EventBroker().register(PeakParamsSuggestedEvent, suggested.append)

    model.suggest_peak_params(make_suggest_request(shape="Pseudo-Voigt"))

    assert len(errors) == 1
    assert suggested == []


def test_suggest_peak_params_non_numeric_range_reports_error(model):
    errors = []
    EventBroker().register(ExceptionEvent, errors.append)

    model.suggest_peak_params(make_suggest_request(range_min="abc"))

    assert len(errors) == 1


def test_suggest_peak_params_range_with_too_few_points_reports_error(model):
    errors = []
    EventBroker().register(ExceptionEvent, errors.append)

    model.suggest_peak_params(make_suggest_request(range_min="4.99", range_max="5.0"))

    assert len(errors) == 1


# ---------------------------------------------------------------------------
# suggest_background_params
# ---------------------------------------------------------------------------

BACKGROUND_SLOPE = 0.4
BACKGROUND_INTERCEPT = 2.0


def make_sloped_background_xy(n=101):
    """A pure straight line, so the guess has an exact answer to recover."""
    x = np.linspace(-5, 5, n)
    y = BACKGROUND_SLOPE * x + BACKGROUND_INTERCEPT
    return x, y


def make_suggest_background_request(**overrides):
    x, y = make_sloped_background_xy()
    defaults = dict(
        source_scan_uuid=UUID(value="scan-001"),
        x=x.tolist(),
        y=y.tolist(),
        range_min="-5",
        range_max="5",
        background="Linear",
    )
    defaults.update(overrides)
    return SuggestBackgroundParamsRequest(**defaults)


def test_suggest_background_params_returns_ok(model):
    response = model.suggest_background_params(make_suggest_background_request())
    assert response.code == ResponseCode.OK


def test_suggest_background_params_publishes_event(model):
    received = []
    EventBroker().register(BackgroundParamsSuggestedEvent, received.append)

    model.suggest_background_params(make_suggest_background_request())

    assert len(received) == 1
    assert received[0].source_scan_uuid == UUID(value="scan-001")


def test_suggest_background_params_recovers_a_straight_line(model):
    received = []
    EventBroker().register(BackgroundParamsSuggestedEvent, received.append)

    model.suggest_background_params(make_suggest_background_request())

    assert received[0].slope == pytest.approx(BACKGROUND_SLOPE, abs=1e-6)
    assert received[0].intercept == pytest.approx(BACKGROUND_INTERCEPT, abs=1e-6)


def test_suggest_background_params_honours_the_fitting_range(model):
    received = []
    EventBroker().register(BackgroundParamsSuggestedEvent, received.append)
    # Outside the range the line bends; a range-respecting guess never sees those points.
    x, y = make_sloped_background_xy()
    y[x > 0] = 100.0

    model.suggest_background_params(make_suggest_background_request(y=y.tolist(), range_min="-5", range_max="0"))

    assert received[0].slope == pytest.approx(BACKGROUND_SLOPE, abs=1e-6)


def test_suggest_background_params_unsupported_background_reports_error(model):
    errors = []
    EventBroker().register(ExceptionEvent, errors.append)
    suggested = []
    EventBroker().register(BackgroundParamsSuggestedEvent, suggested.append)

    model.suggest_background_params(make_suggest_background_request(background="Quadratic"))

    assert len(errors) == 1
    assert suggested == []


def test_suggest_background_params_non_numeric_range_reports_error(model):
    errors = []
    EventBroker().register(ExceptionEvent, errors.append)

    model.suggest_background_params(make_suggest_background_request(range_min="abc"))

    assert len(errors) == 1


def test_perform_fit_publishes_save_fit_event_before_sync_fit_event(model):
    received = []
    EventBroker().register(SaveFitEvent, received.append)
    EventBroker().register(SyncFitEvent, received.append)

    model.perform_fit(make_request())

    assert [type(e) for e in received] == [SaveFitEvent, SyncFitEvent]
    assert received[0].fit_uuid == received[1].fit_uuid
    assert received[0].members == [outcome.member for outcome in received[1].outcomes]


def test_perform_fit_with_every_member_failing_saves_nothing(model):
    saved = []
    EventBroker().register(SaveFitEvent, saved.append)

    model.perform_fit(make_request(range_min="abc"))

    assert saved == []


def test_fit_recompute_syncs_without_saving(model, raw_scans):
    series = make_sequential_scans(raw_scans, count=2)
    entry = make_fit_entry(members=[make_member(s) for s in series])
    saved, synced = [], []
    EventBroker().register(SaveFitEvent, saved.append)
    EventBroker().register(SyncFitEvent, synced.append)

    EventBroker().publish(FitRecomputeEvent(fits=[entry]))

    assert saved == []
    assert len(synced) == 1


# ---------------------------------------------------------------------------
# sync_fit_spec
# ---------------------------------------------------------------------------


def test_sync_fit_spec_publishes_the_saved_member(model, fits):
    member = make_member()
    entry = make_fit_entry(members=[member])
    fits[entry.uuid] = entry
    received = []
    EventBroker().register(SyncFitSpecEvent, received.append)

    response = model.sync_fit_spec(entry.uuid, member.source_scan_uuid)

    assert response.code == ResponseCode.OK
    assert len(received) == 1
    assert received[0].fit_uuid == entry.uuid
    assert received[0].member == member


def test_sync_fit_spec_unknown_fit_is_noop(model):
    received = []
    EventBroker().register(SyncFitSpecEvent, received.append)

    model.sync_fit_spec(UUID(value="fit-999"), UUID(value="scan-001"))

    assert received == []


def test_sync_fit_spec_unknown_member_is_noop(model, fits):
    entry = make_fit_entry()
    fits[entry.uuid] = entry
    received = []
    EventBroker().register(SyncFitSpecEvent, received.append)

    model.sync_fit_spec(entry.uuid, UUID(value="scan-999"))

    assert received == []


def test_sync_fit_spec_reads_the_live_fits_handle(model, fits):
    """A fit saved after the model was built must still be found - the handle is TaviData's own dict."""
    member = make_member()
    entry = make_fit_entry(members=[member])
    received = []
    EventBroker().register(SyncFitSpecEvent, received.append)

    fits[entry.uuid] = entry
    model.sync_fit_spec(entry.uuid, member.source_scan_uuid)

    assert len(received) == 1
