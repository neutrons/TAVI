"""Tests for FitModel."""

import math

import numpy as np
import pytest

from tavi.backend.model.fit_model import FitModel
from tavi.library.data.fit_entry import (
    FitEntry,
    FitRequest,
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
    FitComputedEvent,
    FitRecomputeEvent,
    PeakParamsSuggestedEvent,
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


def make_gaussian_xy(n=101):
    x = np.linspace(-5, 5, n)
    y = GAUSSIAN_AMPLITUDE * np.exp(-((x - GAUSSIAN_CENTER) ** 2) / (2 * GAUSSIAN_SIGMA**2))
    return x, y


def make_series(uuid_val="scan-001") -> PlotSeries:
    return PlotSeries(
        source_scan_uuid=UUID(value=uuid_val),
        scan_name="test_scan",
        normalized_by=None,
        x_name="qh",
        y_name="detector",
        error_name="err",
    )


def make_request(**overrides) -> FitRequest:
    x, y = make_gaussian_xy()
    defaults = dict(
        series=make_series(),
        x=x.tolist(),
        y=y.tolist(),
        err=np.ones_like(x).tolist(),
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
    return FitRequest(**defaults)


@pytest.fixture
def model():
    return FitModel(raw_scans={})


def test_perform_fit_returns_ok(model):
    response = model.perform_fit(make_request())
    assert response.code == ResponseCode.OK


def test_perform_fit_publishes_fit_computed_event(model):
    received = []
    EventBroker().register(FitComputedEvent, received.append)

    model.perform_fit(make_request())

    assert len(received) == 1


def test_perform_fit_does_not_cache_the_curve_on_the_spec(model):
    """The cached FitEntry must not carry the fitted (x, best_fit) points - only its own uuid/series/spec."""
    received = []
    EventBroker().register(FitComputedEvent, received.append)

    model.perform_fit(make_request())

    fit = received[0].fit
    assert not hasattr(fit, "x")
    assert not hasattr(fit, "best_fit")
    assert fit.series.source_scan_uuid == UUID(value="scan-001")


def test_perform_fit_curve_carries_its_source_scan(model):
    """The curve identifies the series it belongs to, so redrawing it replaces that series' line."""
    received = []
    EventBroker().register(FitComputedEvent, received.append)

    model.perform_fit(make_request())

    assert received[0].curve.source_scan_uuid == UUID(value="scan-001")


def test_perform_fit_without_a_fit_uuid_mints_a_new_one_each_time(model):
    received = []
    EventBroker().register(FitComputedEvent, received.append)

    model.perform_fit(make_request())
    model.perform_fit(make_request())

    assert received[0].fit.uuid != received[1].fit.uuid


def test_perform_fit_with_a_fit_uuid_refits_that_fit_in_place(model):
    """A refit carries the existing uuid through, so TaviProjectModel overwrites that entry."""
    received = []
    EventBroker().register(FitComputedEvent, received.append)

    model.perform_fit(make_request(fit_uuid=UUID(value="fit-007")))

    assert received[0].fit.uuid == UUID(value="fit-007")


def test_perform_fit_recovers_amplitude(model):
    received = []
    EventBroker().register(FitComputedEvent, received.append)

    model.perform_fit(make_request())

    curve = received[0].curve
    assert max(curve.best_fit) == pytest.approx(GAUSSIAN_AMPLITUDE, rel=1e-3)


def test_perform_fit_reduced_chi_squared_near_zero_for_noiseless_data(model):
    received = []
    EventBroker().register(FitComputedEvent, received.append)

    model.perform_fit(make_request())

    assert received[0].result.reduced_chi_squared == pytest.approx(0.0, abs=1e-6)


def test_perform_fit_result_recovers_center_and_fwhm(model):
    received = []
    EventBroker().register(FitComputedEvent, received.append)

    model.perform_fit(make_request())

    result = received[0].result
    assert result.amplitude == pytest.approx(GAUSSIAN_AREA, rel=1e-3)
    assert result.center == pytest.approx(GAUSSIAN_CENTER, abs=1e-3)
    assert result.fwhm == pytest.approx(GAUSSIAN_FWHM, rel=1e-3)


def test_perform_fit_reports_uncertainties_for_noisy_data(model):
    """Perfectly noiseless synthetic data can leave lmfit's covariance degenerate (stderr=None) -
    add a touch of deterministic noise so there's something real to estimate an uncertainty from."""
    received = []
    EventBroker().register(FitComputedEvent, received.append)

    x, y = make_gaussian_xy()
    rng = np.random.default_rng(seed=0)
    noisy_y = (y + rng.normal(scale=0.02, size=y.shape)).tolist()
    model.perform_fit(make_request(y=noisy_y, err=np.full_like(x, 0.02).tolist()))

    result = received[0].result
    assert result.amplitude_err is not None
    assert result.center_err is not None
    assert result.fwhm_err is not None


def test_perform_fit_with_constant_background_succeeds(model):
    received = []
    EventBroker().register(FitComputedEvent, received.append)

    x, y = make_gaussian_xy()
    request = make_request(background="Constant", background_constant=make_param(0), y=(y + 2.0).tolist())

    model.perform_fit(request)

    assert len(received) == 1
    assert received[0].result.background_constant == pytest.approx(2.0, abs=1e-2)


def test_perform_fit_trims_to_range(model):
    received = []
    EventBroker().register(FitComputedEvent, received.append)

    model.perform_fit(make_request(range_min="-1", range_max="1"))

    curve = received[0].curve
    assert min(curve.x) >= -1
    assert max(curve.x) <= 1
    assert len(curve.x) < 101


def test_perform_fit_unsupported_peak_shape_reports_error(model):
    errors = []
    EventBroker().register(ExceptionEvent, errors.append)
    computed = []
    EventBroker().register(FitComputedEvent, computed.append)

    request = make_request(peak=PeakField(shape="Pseudo-Voigt", amplitude=make_param(1), center=make_param(0), fwhm=make_param(1)))
    model.perform_fit(request)

    assert len(errors) == 1
    assert computed == []


def test_perform_fit_unsupported_background_reports_error(model):
    errors = []
    EventBroker().register(ExceptionEvent, errors.append)
    computed = []
    EventBroker().register(FitComputedEvent, computed.append)

    model.perform_fit(make_request(background="Quadratic"))

    assert len(errors) == 1
    assert computed == []


def test_perform_fit_non_numeric_peak_field_reports_error(model):
    errors = []
    EventBroker().register(ExceptionEvent, errors.append)

    request = make_request(
        peak=PeakField(shape="Gaussian", amplitude=make_param("not-a-number"), center=make_param(0), fwhm=make_param(1))
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


# ---------------------------------------------------------------------------
# Blank initial params -> auto-guess from data
# ---------------------------------------------------------------------------


def test_perform_fit_guesses_blank_peak_params(model):
    """Leaving all three peak fields blank must guess from the data, not error out."""
    received = []
    EventBroker().register(FitComputedEvent, received.append)
    errors = []
    EventBroker().register(ExceptionEvent, errors.append)

    request = make_request(
        peak=PeakField(shape="Gaussian", amplitude=make_param(""), center=make_param(""), fwhm=make_param(""))
    )
    model.perform_fit(request)

    assert errors == []
    assert len(received) == 1
    assert received[0].result.amplitude == pytest.approx(GAUSSIAN_AREA, rel=1e-2)
    assert received[0].result.center == pytest.approx(GAUSSIAN_CENTER, abs=1e-2)


def test_perform_fit_guesses_only_the_blank_peak_field(model):
    """A deliberately-wrong explicit center must survive - only the blank fields get guessed."""
    received = []
    EventBroker().register(FitComputedEvent, received.append)

    request = make_request(
        peak=PeakField(
            shape="Gaussian", amplitude=make_param(""), center=make_param(GAUSSIAN_CENTER), fwhm=make_param("")
        )
    )
    model.perform_fit(request)

    assert received[0].result.center == pytest.approx(GAUSSIAN_CENTER, abs=1e-3)


def test_perform_fit_guesses_blank_background_constant(model):
    received = []
    EventBroker().register(FitComputedEvent, received.append)

    x, y = make_gaussian_xy()
    request = make_request(background="Constant", background_constant=make_param(""), y=(y + 2.0).tolist())
    model.perform_fit(request)

    assert received[0].result.background_constant == pytest.approx(2.0, abs=0.5)


# ---------------------------------------------------------------------------
# FitRecomputeEvent - recompute a selected fit against live data, never a cached curve
# ---------------------------------------------------------------------------


def make_fit_entry_from_request(request: FitRequest, uuid_val="fit-001"):
    return FitEntry(
        uuid=UUID(value=uuid_val),
        series=request.series,
        range_min=request.range_min,
        range_max=request.range_max,
        background=request.background,
        background_constant=request.background_constant,
        peak=request.peak,
    )


def make_scan(uuid_val="scan-001"):
    x, y = make_gaussian_xy()
    return RawScan(
        uuid=UUID(value=uuid_val),
        data=ScanData(data={"qh": x.tolist(), "detector": y.tolist()}),
        metadata=ScanMetadata(),
        tavimeta=TaviMetadata(default_axis=("qh", "detector"), friendly_name="test_scan", friendly_path="/exp1"),
        prov=Provenance(raw_file="scan.dat", contributing_scans={UUID(value=uuid_val): 1}),
    )


def test_fit_focus_recomputes_against_current_raw_scan_data():
    scan = make_scan()
    raw_scans = {scan.uuid: scan}
    FitModel(raw_scans=raw_scans)
    entry = make_fit_entry_from_request(make_request())

    received = []
    EventBroker().register(FitComputedEvent, received.append)
    EventBroker().publish(FitRecomputeEvent(fits=[entry]))

    assert len(received) == 1
    assert received[0].fit.uuid == entry.uuid
    assert received[0].result.amplitude == pytest.approx(GAUSSIAN_AREA, rel=1e-2)


def test_fit_focus_reflects_changed_underlying_data():
    """Recompute must reflect the scan's *current* data, not whatever it was fit against originally."""
    scan = make_scan()
    raw_scans = {scan.uuid: scan}
    FitModel(raw_scans=raw_scans)
    entry = make_fit_entry_from_request(make_request())

    doubled_amplitude = (2 * GAUSSIAN_AMPLITUDE * np.exp(-(np.array(scan.data.data["qh"]) ** 2) / 2)).tolist()
    raw_scans[scan.uuid].data.data["detector"] = doubled_amplitude

    received = []
    EventBroker().register(FitComputedEvent, received.append)
    EventBroker().publish(FitRecomputeEvent(fits=[entry]))

    assert received[0].result.amplitude == pytest.approx(2 * GAUSSIAN_AREA, rel=1e-2)


def test_fit_focus_unknown_source_scan_is_noop():
    FitModel(raw_scans={})
    entry = make_fit_entry_from_request(make_request())

    received = []
    EventBroker().register(FitComputedEvent, received.append)
    EventBroker().publish(FitRecomputeEvent(fits=[entry]))

    assert received == []


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


def make_linear_background_request(**overrides) -> FitRequest:
    x, y = make_gaussian_on_a_slope_xy()
    defaults = dict(
        x=x.tolist(),
        y=y.tolist(),
        err=np.ones_like(x).tolist(),
        background="Linear",
        background_constant=make_param(""),
        background_slope=make_param(""),
    )
    defaults.update(overrides)
    return make_request(**defaults)


def test_perform_fit_linear_background_is_supported(model):
    errors = []
    EventBroker().register(ExceptionEvent, errors.append)
    computed = []
    EventBroker().register(FitComputedEvent, computed.append)

    model.perform_fit(make_linear_background_request())

    assert errors == []
    assert len(computed) == 1


def test_perform_fit_linear_background_recovers_slope_and_intercept(model):
    computed = []
    EventBroker().register(FitComputedEvent, computed.append)

    model.perform_fit(make_linear_background_request())

    result = computed[0].result
    assert result.background_slope == pytest.approx(LINEAR_BG_SLOPE, abs=0.05)
    assert result.background_constant == pytest.approx(LINEAR_BG_INTERCEPT, abs=0.1)


def test_perform_fit_linear_background_still_recovers_the_peak(model):
    computed = []
    EventBroker().register(FitComputedEvent, computed.append)

    model.perform_fit(make_linear_background_request())

    result = computed[0].result
    assert result.center == pytest.approx(GAUSSIAN_CENTER, abs=0.05)
    assert result.amplitude == pytest.approx(GAUSSIAN_AREA, rel=0.05)


def test_perform_fit_constant_background_reports_no_slope(model):
    computed = []
    EventBroker().register(FitComputedEvent, computed.append)

    model.perform_fit(make_request(background="Constant"))

    assert computed[0].result.background_slope is None


def test_perform_fit_linear_background_honours_a_fixed_slope(model):
    computed = []
    EventBroker().register(FitComputedEvent, computed.append)

    model.perform_fit(make_linear_background_request(background_slope=make_param(0, fixed=True)))

    assert computed[0].result.background_slope == pytest.approx(0.0, abs=1e-9)


def test_perform_fit_linear_background_non_numeric_slope_reports_error(model):
    errors = []
    EventBroker().register(ExceptionEvent, errors.append)
    computed = []
    EventBroker().register(FitComputedEvent, computed.append)

    model.perform_fit(make_linear_background_request(background_slope=make_param("abc")))

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


def test_suggest_background_params_supports_constant(model):
    received = []
    EventBroker().register(BackgroundParamsSuggestedEvent, received.append)

    model.suggest_background_params(make_suggest_background_request(background="Constant"))

    assert len(received) == 1


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
