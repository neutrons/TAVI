"""Tests for GroupFitModel."""

import math

import numpy as np
import pytest

from tavi.backend.model.fit_model import FitModel
from tavi.backend.model.group_fit_model import GroupFitModel, _group_fit_params, _PreparedScan
from tavi.library.data.fit_entry import CosineField, FitResultSummary, FitSpec, ParamField, PeakField, PeakResult
from tavi.library.data.group_fit_entry import (
    GroupFitRequest,
    SetGroupTitleRequest,
    SuggestGroupBackgroundParamsRequest,
    SuggestGroupPeakParamsRequest,
    SummarizeGroupFitRequest,
)
from tavi.library.data.model_response import ResponseCode
from tavi.library.data.scan import UUID, ProcessedScan, Provenance, RawScan, ScanData, ScanMetadata, TaviMetadata
from tavi.library.fit.fit import ModelName
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.exception_event import ReportErrorEvent
from tavi.meta.event.type.presenter_event import (
    RecomputeFitEvent,
    SaveFitSummaryEvent,
    SyncBackgroundParamsEvent,
    SyncCosineParamsEvent,
    SyncGroupBackgroundParamsEvent,
    SyncGroupCosineParamsEvent,
    SyncGroupFitEvent,
    SyncGroupPeakParamsEvent,
    SyncGroupTitlesEvent,
    SyncPeakParamsEvent,
)

AMPLITUDE = 50.0
SIGMA = 0.5
FWHM = SIGMA * 2 * math.sqrt(2 * math.log(2))
# A temperature sweep: the peak drifts one step further along x in each scan.
CENTERS = [-1.0, 0.0, 1.0]


def make_param(value="", fixed=False, minimum="", maximum="", constrained=False) -> ParamField:
    return ParamField(
        value=str(value), fixed=fixed, minimum=str(minimum), maximum=str(maximum), constrained=constrained
    )


def make_peak(shape="Gaussian", amplitude="", center="", fwhm="") -> PeakField:
    return PeakField(shape=shape, amplitude=make_param(amplitude), center=make_param(center), fwhm=make_param(fwhm))


def make_spec(peaks=None, **overrides) -> FitSpec:
    fields = dict(
        range_min="-5",
        range_max="5",
        background="None",
        background_constant=make_param(0),
        peaks=peaks if peaks is not None else [make_peak(amplitude=AMPLITUDE, center=CENTERS[0], fwhm=FWHM)],
    )
    fields.update(overrides)
    return FitSpec(**fields)


def make_columns(center, n=201, temp=10.0):
    x = np.linspace(-5, 5, n)
    detector = AMPLITUDE * np.exp(-((x - center) ** 2) / (2 * SIGMA**2)) + 1.0
    return {"qh": x.tolist(), "detector": detector.tolist(), "temp": np.full(n, temp).tolist()}


def make_metadata(name):
    return TaviMetadata(default_axis=("qh", "detector"), friendly_name=name, friendly_path="/exp1")


def make_scan(uuid_val, center=0.0, temp=10.0, columns=None) -> RawScan:
    return RawScan(
        uuid=UUID(value=uuid_val),
        data=ScanData(data=columns if columns is not None else make_columns(center, temp=temp)),
        metadata=ScanMetadata(),
        tavimeta=make_metadata(uuid_val),
        prov=Provenance(raw_file=f"{uuid_val}.dat", contributing_scans={UUID(value=uuid_val): 1}),
    )


def make_processed_scan(uuid_val, center=0.0) -> ProcessedScan:
    return ProcessedScan(
        uuid=UUID(value=uuid_val),
        data=ScanData(data=make_columns(center)),
        metadata=ScanMetadata(),
        tavimeta=make_metadata(uuid_val),
        prov=Provenance(raw_file="", contributing_scans={UUID(value=uuid_val): 1}),
    )


def make_sweep(temps=(10.0, 20.0, 30.0)) -> dict:
    return {
        UUID(value=f"s{i}"): make_scan(f"s{i}", center=center, temp=temp)
        for i, (center, temp) in enumerate(zip(CENTERS, temps))
    }


def make_request(specs=None, scans=None, **overrides) -> GroupFitRequest:
    uuids = list(scans) if scans is not None else [UUID(value=f"s{i}") for i in range(len(CENTERS))]
    fields = dict(specs=specs if specs is not None else {uuid: make_spec() for uuid in uuids}, y_col="detector")
    fields.update(overrides)
    return GroupFitRequest(**fields)


def make_prepared(uuid_val, spec) -> _PreparedScan:
    """A scan prepared the way perform_group_fit does, for testing _group_fit_params on its own."""
    scan = make_scan(uuid_val)
    x = np.asarray(scan.data.data["qh"])
    y = np.asarray(scan.data.data["detector"])
    model_dict = GroupFitModel({}, {})._translator._build_model_dict(spec, x, y)
    return _PreparedScan(uuid=scan.uuid, name=uuid_val, spec=spec, x=x, y=y, model_dict=model_dict, trimmed_scan=scan)


@pytest.fixture
def raw_scans():
    return make_sweep()


@pytest.fixture
def model(raw_scans):
    return GroupFitModel(raw_scans, {})


@pytest.fixture
def synced():
    received = []
    EventBroker().register(SyncGroupFitEvent, received.append)
    return received


@pytest.fixture
def reported():
    received = []
    EventBroker().register(ReportErrorEvent, received.append)
    return received


def centers(event):
    return [outcome.result.peaks[0].center for outcome in event.outcomes]


# ---------------------------------------------------------------------------
# Wiring
# ---------------------------------------------------------------------------


def test_group_fit_model_does_not_subscribe_to_recompute():
    """Borrowing FitModel's helpers must not add a second recompute of every saved fit."""
    before = list(EventBroker().registry.get(RecomputeFitEvent, []))

    GroupFitModel({}, {})

    assert EventBroker().registry.get(RecomputeFitEvent, []) == before


def test_fit_model_still_subscribes_to_recompute():
    fit_model = FitModel({}, {})

    assert fit_model._handle_fit_recompute_event in EventBroker().registry[RecomputeFitEvent]


# ---------------------------------------------------------------------------
# perform_group_fit
# ---------------------------------------------------------------------------


def test_publishes_one_sync_with_an_outcome_per_scan_in_order(model, synced):
    response = model.perform_group_fit(make_request())

    assert response.code == ResponseCode.OK
    assert len(synced) == 1
    assert [outcome.source_scan_uuid for outcome in synced[0].outcomes] == [UUID(value=f"s{i}") for i in range(3)]
    assert (synced[0].x_col, synced[0].y_col) == (None, "detector")


def test_chained_fit_tracks_the_drifting_peak(model, synced):
    model.perform_group_fit(make_request(chain=True))

    assert centers(synced[0]) == pytest.approx(CENTERS, abs=1e-3)


def test_unchained_fit_from_each_scans_own_spec(model, synced):
    specs = {UUID(value=f"s{i}"): make_spec(peaks=[make_peak(center=c)]) for i, c in enumerate(CENTERS)}

    model.perform_group_fit(make_request(specs=specs, chain=False))

    assert centers(synced[0]) == pytest.approx(CENTERS, abs=1e-3)


def test_outcomes_carry_spec_curve_and_data(model, synced):
    request = make_request()

    model.perform_group_fit(request)

    outcome = synced[0].outcomes[0]
    assert outcome.spec == request.specs[outcome.source_scan_uuid]
    assert len(outcome.curve.x) == 300
    assert outcome.curve.scan_name == "s0"
    assert outcome.data.err == pytest.approx(np.sqrt(np.abs(outcome.data.y)).tolist())
    assert outcome.result.reduced_chi_squared >= 0


def test_each_scans_range_trims_its_own_data(model, synced):
    specs = {uuid: make_spec() for uuid in make_sweep()}
    specs[UUID(value="s1")] = make_spec(range_min="-2", range_max="2")

    model.perform_group_fit(make_request(specs=specs))

    trimmed = synced[0].outcomes[1].data.x
    assert min(trimmed) >= -2 and max(trimmed) <= 2
    assert min(synced[0].outcomes[0].data.x) == pytest.approx(-5)


def test_averages_cover_the_whole_scan_not_its_range(raw_scans, synced):
    columns = make_columns(CENTERS[0])
    # A temperature that changes only outside the fitting range must still count towards its mean.
    columns["temp"] = [100.0 if x < -2 else 10.0 for x in columns["qh"]]
    raw_scans[UUID(value="s0")] = make_scan("s0", columns=columns)
    specs = {uuid: make_spec(range_min="-2", range_max="5") for uuid in raw_scans}

    GroupFitModel(raw_scans, {}).perform_group_fit(make_request(specs=specs, averages=["temp"]))

    averages = [outcome.averages["temp"] for outcome in synced[0].outcomes]
    assert averages[0].mean == pytest.approx(np.mean(columns["temp"]))
    assert [a.mean for a in averages[1:]] == pytest.approx([20.0, 30.0])
    assert averages[1].std == pytest.approx(0.0)


def test_explicit_x_col_is_used(raw_scans, synced):
    for uuid, scan in list(raw_scans.items()):
        columns = {**scan.data.data, "q2": scan.data.data["qh"]}
        raw_scans[uuid] = make_scan(uuid.value, columns=columns)
    model = GroupFitModel(raw_scans, {})

    model.perform_group_fit(make_request(x_col="q2"))

    assert synced[0].x_col == "q2"
    assert centers(synced[0]) == pytest.approx(CENTERS, abs=1e-3)


def test_digit_leading_columns_find_their_underscore_spelling(raw_scans, synced):
    for uuid, scan in list(raw_scans.items()):
        columns = {**scan.data.data, "_2theta": scan.data.data["qh"]}
        raw_scans[uuid] = make_scan(uuid.value, columns=columns)

    GroupFitModel(raw_scans, {}).perform_group_fit(make_request(x_col="2theta"))

    assert centers(synced[0]) == pytest.approx(CENTERS, abs=1e-3)


def test_processed_scans_fit_too(synced):
    processed = {UUID(value="c0"): make_processed_scan("c0", center=0.5)}

    GroupFitModel({}, processed).perform_group_fit(
        make_request(specs={UUID(value="c0"): make_spec(peaks=[make_peak(center=0.5)])})
    )

    assert centers(synced[0]) == pytest.approx([0.5], abs=1e-3)


def test_a_linear_background_is_fit(model, synced):
    spec = make_spec(background="Linear", background_constant=make_param(""), background_slope=make_param(""))

    model.perform_group_fit(make_request(specs={uuid: spec for uuid in make_sweep()}))

    assert [outcome.result.background_constant for outcome in synced[0].outcomes] == pytest.approx(
        [1.0, 1.0, 1.0], abs=1e-3
    )


def test_does_not_modify_the_project_scans(model, raw_scans):
    before = {uuid: list(scan.data.data["qh"]) for uuid, scan in raw_scans.items()}
    specs = {uuid: make_spec(range_min="-1", range_max="1") for uuid in raw_scans}

    model.perform_group_fit(make_request(specs=specs))

    assert {uuid: list(scan.data.data["qh"]) for uuid, scan in raw_scans.items()} == before


# ---------------------------------------------------------------------------
# Scans left out, and failures
# ---------------------------------------------------------------------------


def test_an_unloaded_scan_is_reported_and_listed_without_result(model, synced, reported):
    specs = {uuid: make_spec() for uuid in make_sweep()}
    specs[UUID(value="gone")] = make_spec()

    model.perform_group_fit(make_request(specs=specs))

    assert len(reported) == 1
    assert "gone: Scan is no longer loaded." in reported[0].error.message
    outcomes = synced[0].outcomes
    assert outcomes[-1].source_scan_uuid == UUID(value="gone") and outcomes[-1].result is None
    assert all(outcome.result is not None for outcome in outcomes[:-1])


def test_a_scan_whose_model_differs_is_left_out(model, synced, reported):
    specs = {uuid: make_spec() for uuid in make_sweep()}
    specs[UUID(value="s1")] = make_spec(peaks=[make_peak(), make_peak(center=2)])

    model.perform_group_fit(make_request(specs=specs))

    assert len(reported) == 1
    assert "s1: Its background and peaks differ from s0's" in reported[0].error.message
    assert [outcome.result is None for outcome in synced[0].outcomes] == [False, True, False]


def test_a_bad_range_is_reported_and_left_out(model, synced, reported):
    specs = {uuid: make_spec() for uuid in make_sweep()}
    specs[UUID(value="s2")] = make_spec(range_min="low")

    model.perform_group_fit(make_request(specs=specs))

    assert "s2: Fitting range 'low'" in reported[0].error.message
    assert synced[0].outcomes[2].result is None


def test_problems_are_gathered_into_one_report(model, reported):
    specs = {uuid: make_spec(range_min="low") for uuid in make_sweep()}

    model.perform_group_fit(make_request(specs=specs))

    assert len(reported) == 1
    assert reported[0].error.message.count("Fitting range") == 3


def test_a_missing_y_column_publishes_nothing(model, synced, reported):
    model.perform_group_fit(make_request(y_col="monitor"))

    assert synced == []
    assert "Column 'monitor' not found" in reported[0].error.message


def test_a_missing_average_column_publishes_nothing(model, synced, reported):
    model.perform_group_fit(make_request(averages=["field"]))

    assert synced == []
    assert "Group fit: Column 'field' not found" in reported[0].error.message


def test_an_empty_request_is_reported(model, synced, reported):
    model.perform_group_fit(make_request(specs={}))

    assert synced == []
    assert "A group fit needs at least one scan." in reported[0].error.message


# ---------------------------------------------------------------------------
# _group_fit_params - request -> GroupFitParams
# ---------------------------------------------------------------------------


def test_params_use_the_first_scans_model_for_the_group():
    first, second = make_prepared("s0", make_spec()), make_prepared("s1", make_spec())

    params = _group_fit_params(make_request(x_col="qh", averages=["temp"]), [first, second])

    assert params.fits == {"detector": first.model_dict}
    assert (params.x_col, params.averages, params.chain) == ("qh", ["temp"], True)


def test_chained_params_override_nothing_a_scan_did_not_change():
    params = _group_fit_params(make_request(), [make_prepared("s0", make_spec()), make_prepared("s1", make_spec())])

    assert params.overrides == {}


def test_chained_params_override_only_the_changed_parameter():
    edited = make_spec(peaks=[make_peak(amplitude=AMPLITUDE, center=0.25, fwhm=FWHM)])
    edited.peaks[0].center.fixed = True

    params = _group_fit_params(make_request(), [make_prepared("s0", make_spec()), make_prepared("s1", edited)])

    assert params.overrides == {
        UUID(value="s1"): {"detector": {"peak1_": {"center": 0.25, "set": {"center": {"vary": False}}}}}
    }


def test_chained_params_override_fwhm_as_lmfit_sigma():
    edited = make_spec(peaks=[make_peak(amplitude=AMPLITUDE, center=CENTERS[0], fwhm=2 * FWHM)])

    params = _group_fit_params(make_request(), [make_prepared("s0", make_spec()), make_prepared("s1", edited)])

    patch = params.overrides[UUID(value="s1")]["detector"]["peak1_"]
    assert patch["sigma"] == pytest.approx(2 * SIGMA)


def test_chained_params_leave_a_blank_changed_value_to_the_seed():
    """A field emptied but bounded still takes its start from the previous scan, only its bounds change."""
    edited = make_spec()
    edited.peaks[0].center = make_param("", minimum="-3", maximum="3")

    params = _group_fit_params(make_request(), [make_prepared("s0", make_spec()), make_prepared("s1", edited)])

    patch = params.overrides[UUID(value="s1")]["detector"]["peak1_"]
    assert patch == {"set": {"center": {"vary": True, "min": -3.0, "max": 3.0}}}


def test_chained_params_compare_cosine_parameters():
    def cosine_spec(frequency):
        fields = {name: make_param(1) for name in ("amplitude", "phase", "offset")}
        return make_spec(
            peaks=[PeakField(shape="Cosine", cosine=CosineField(frequency=make_param(frequency), **fields))]
        )

    params = _group_fit_params(
        make_request(), [make_prepared("s0", cosine_spec(1)), make_prepared("s1", cosine_spec(2))]
    )

    assert params.overrides[UUID(value="s1")]["detector"] == {
        "peak1_": {"frequency": 2.0, "set": {"frequency": {"vary": True}}}
    }


def test_chained_params_compare_background_parameters():
    def linear(slope):
        return make_spec(background="Linear", background_slope=make_param(slope))

    params = _group_fit_params(make_request(), [make_prepared("s0", linear(0)), make_prepared("s1", linear(0.5))])

    assert params.overrides[UUID(value="s1")]["detector"] == {"bg_": {"slope": 0.5, "set": {"slope": {"vary": True}}}}


def test_unchained_params_override_every_component_of_every_later_scan():
    first, second = make_prepared("s0", make_spec()), make_prepared("s1", make_spec())

    params = _group_fit_params(make_request(chain=False), [first, second])

    ((model_name, component),) = second.model_dict
    assert model_name is ModelName.Gaussian
    assert params.overrides == {
        UUID(value="s1"): {"detector": {"peak1_": {k: v for k, v in component.items() if k != "prefix"}}}
    }
    assert params.chain is False


def test_a_blank_y_column_is_reported(model, synced, reported):
    model.perform_group_fit(make_request(y_col="  "))

    assert synced == []
    assert reported[0].error.message == "Set the y column to fit."


def test_a_blank_range_fits_the_whole_scan(model, synced):
    specs = {uuid: make_spec(range_min="", range_max=" ") for uuid in make_sweep()}

    model.perform_group_fit(make_request(specs=specs))

    data_x = synced[0].outcomes[0].data.x
    assert (min(data_x), max(data_x), len(data_x)) == (pytest.approx(-5), pytest.approx(5), 201)


def test_a_blank_bound_leaves_only_that_side_open(model, synced):
    specs = {uuid: make_spec(range_min="", range_max="2") for uuid in make_sweep()}

    model.perform_group_fit(make_request(specs=specs))

    data_x = synced[0].outcomes[0].data.x
    assert min(data_x) == pytest.approx(-5) and max(data_x) <= 2


# ---------------------------------------------------------------------------
# Suggest Params.
# ---------------------------------------------------------------------------


def make_peak_guess(uuid_val="s0", shape="Gaussian", **overrides) -> SuggestGroupPeakParamsRequest:
    fields = dict(source_scan_uuid=UUID(value=uuid_val), y_col="detector", range_min="", range_max="", shape=shape)
    fields.update(overrides)
    return SuggestGroupPeakParamsRequest(**fields)


def make_background_guess(uuid_val="s0", **overrides) -> SuggestGroupBackgroundParamsRequest:
    fields = dict(
        source_scan_uuid=UUID(value=uuid_val), y_col="detector", range_min="", range_max="", background="Linear"
    )
    fields.update(overrides)
    return SuggestGroupBackgroundParamsRequest(**fields)


def collect(event_type):
    received = []
    EventBroker().register(event_type, received.append)
    return received


@pytest.mark.parametrize(("uuid_val", "center"), [("s0", CENTERS[0]), ("s2", CENTERS[2])])
def test_peak_guess_is_made_from_that_scans_own_columns(model, uuid_val, center):
    guesses = collect(SyncGroupPeakParamsEvent)

    model.suggest_peak_params(make_peak_guess(uuid_val))

    assert len(guesses) == 1
    assert guesses[0].source_scan_uuid == UUID(value=uuid_val)
    assert guesses[0].center == pytest.approx(center, abs=0.1)


def test_group_guesses_never_reach_the_1d_fitting_tab(model):
    """The tab listens for the single-fit events; a guess asked for in the group window must not land there."""
    singles = [collect(t) for t in (SyncPeakParamsEvent, SyncCosineParamsEvent, SyncBackgroundParamsEvent)]

    model.suggest_peak_params(make_peak_guess())
    model.suggest_peak_params(make_peak_guess(shape="Cosine"))
    model.suggest_background_params(make_background_guess())

    assert all(received == [] for received in singles)


def test_a_cosine_guess_is_a_group_cosine_event(model):
    guesses = collect(SyncGroupCosineParamsEvent)

    model.suggest_peak_params(make_peak_guess(shape="Cosine"))

    assert [guess.source_scan_uuid for guess in guesses] == [UUID(value="s0")]


def test_a_background_guess_is_a_group_background_event(model):
    guesses = collect(SyncGroupBackgroundParamsEvent)

    model.suggest_background_params(make_background_guess("s1"))

    assert len(guesses) == 1 and guesses[0].source_scan_uuid == UUID(value="s1")
    assert np.isfinite([guesses[0].slope, guesses[0].intercept]).all()


def test_a_guess_uses_the_scans_range(model):
    guesses = collect(SyncGroupPeakParamsEvent)

    # s0's peak sits at -1; guessed over 0..5 only, the guess can't be there.
    model.suggest_peak_params(make_peak_guess(range_min="0", range_max="5"))

    assert guesses[0].center >= 0


def test_a_guess_uses_an_explicit_x_column(raw_scans):
    for uuid, scan in list(raw_scans.items()):
        raw_scans[uuid] = make_scan(uuid.value, columns={**scan.data.data, "q2": scan.data.data["qh"]})
    guesses = collect(SyncGroupPeakParamsEvent)

    GroupFitModel(raw_scans, {}).suggest_peak_params(make_peak_guess(x_col="q2"))

    assert guesses[0].center == pytest.approx(CENTERS[0], abs=0.1)


@pytest.mark.parametrize(
    ("request_", "message"),
    [
        (make_peak_guess("gone"), "Scan is no longer loaded."),
        (make_peak_guess(y_col="monitor"), "Column 'monitor' not found"),
        (make_peak_guess(y_col=""), "Set the y column to fit."),
        (make_peak_guess(shape="Custom"), "Peak shape 'Custom' is not supported yet."),
        (make_peak_guess(range_min="low"), "Fitting range 'low'"),
    ],
)
def test_a_guess_that_cannot_be_made_is_reported(model, reported, request_, message):
    guesses = collect(SyncGroupPeakParamsEvent)

    model.suggest_peak_params(request_)

    assert guesses == []
    assert message in reported[0].error.message


def test_an_unsupported_background_guess_is_reported(model, reported):
    model.suggest_background_params(make_background_guess(background="None"))

    assert "Background 'None' is not supported yet." in reported[0].error.message


# ---------------------------------------------------------------------------
# summarize_group_fit
# ---------------------------------------------------------------------------


def make_result(center, center_err=0.01) -> FitResultSummary:
    peak = PeakResult(amplitude=AMPLITUDE, amplitude_err=1.0, center=center, center_err=center_err, fwhm=FWHM)
    return FitResultSummary(reduced_chi_squared=1.0, peaks=[peak])


def make_summarize_request(results=None, **overrides) -> SummarizeGroupFitRequest:
    if results is None:
        results = {UUID(value=f"s{i}"): make_result(center) for i, center in enumerate(CENTERS)}
    fields = dict(parameter="peak1 center", x_col="temp", results=results)
    fields.update(overrides)
    return SummarizeGroupFitRequest(**fields)


@pytest.fixture
def saved():
    received = []
    EventBroker().register(SaveFitSummaryEvent, received.append)
    return received


def test_summarize_saves_the_parameter_against_each_scans_average(model, saved, reported):
    response = model.summarize_group_fit(make_summarize_request())

    assert response.code == ResponseCode.OK
    assert reported == []
    assert len(saved) == 1
    summary = saved[0].summary
    assert (summary.name, summary.x_label, summary.y_label) == ("peak1 center vs temp", "temp", "peak1 center")
    assert [point.scan_name for point in summary.points] == ["s0", "s1", "s2"]
    assert [point.x for point in summary.points] == [10.0, 20.0, 30.0]
    assert [point.y for point in summary.points] == CENTERS
    assert [point.y_err for point in summary.points] == [0.01] * 3


def test_summarize_takes_the_columns_spread_as_its_uncertainty(saved):
    columns = make_columns(0.0)
    columns["temp"] = np.linspace(9.0, 11.0, len(columns["qh"])).tolist()
    model = GroupFitModel({UUID(value="s0"): make_scan("s0", columns=columns)}, {})

    model.summarize_group_fit(make_summarize_request(results={UUID(value="s0"): make_result(0.0)}))

    point = saved[0].summary.points[0]
    assert point.x == pytest.approx(10.0)
    assert point.x_err == pytest.approx(np.std(columns["temp"]))


def test_summarize_a_group_fits_own_results(model, synced, saved):
    model.perform_group_fit(make_request())
    results = {outcome.source_scan_uuid: outcome.result for outcome in synced[0].outcomes}

    model.summarize_group_fit(make_summarize_request(results=results))

    assert [point.y for point in saved[0].summary.points] == pytest.approx(CENTERS, abs=1e-3)


def test_summarize_reads_a_processed_scan_too(saved):
    model = GroupFitModel({}, {UUID(value="c0"): make_processed_scan("c0")})

    model.summarize_group_fit(make_summarize_request(results={UUID(value="c0"): make_result(0.0)}))

    assert saved[0].summary.points[0].x == 10.0


def test_summarize_without_a_column_reports_and_saves_nothing(model, saved, reported):
    model.summarize_group_fit(make_summarize_request(x_col="  "))

    assert saved == []
    assert len(reported) == 1


def test_summarize_without_results_reports_and_saves_nothing(model, saved, reported):
    model.summarize_group_fit(make_summarize_request(results={}))

    assert saved == []
    assert len(reported) == 1


def test_summarize_leaves_out_a_scan_missing_the_column(raw_scans, saved, reported):
    columns = make_columns(0.0)
    del columns["temp"]
    raw_scans[UUID(value="s1")] = make_scan("s1", columns=columns)

    GroupFitModel(raw_scans, {}).summarize_group_fit(make_summarize_request())

    assert [point.scan_name for point in saved[0].summary.points] == ["s0", "s2"]
    assert len(reported) == 1
    assert "s1" in reported[0].error.args[0]


def test_summarize_leaves_out_a_scan_no_longer_loaded(model, saved, reported):
    results = {UUID(value="s0"): make_result(0.0), UUID(value="gone"): make_result(1.0)}

    model.summarize_group_fit(make_summarize_request(results=results))

    assert [point.scan_name for point in saved[0].summary.points] == ["s0"]
    assert len(reported) == 1


def test_summarize_an_unknown_parameter_reports_once_and_saves_nothing(model, saved, reported):
    model.summarize_group_fit(make_summarize_request(parameter="peak2 center"))

    assert saved == []
    assert len(reported) == 1


# ---------------------------------------------------------------------------
# Set Title
# ---------------------------------------------------------------------------


def make_titled_scan(uuid_val, **metadata) -> RawScan:
    return make_scan(uuid_val).model_copy(update={"metadata": ScanMetadata(data=metadata)})


def set_title(scans, title, uuids=None) -> list:
    received = []
    EventBroker().register(SyncGroupTitlesEvent, received.append)
    model = GroupFitModel({scan.uuid: scan for scan in scans}, {})
    model.set_title(SetGroupTitleRequest(title=title, scan_uuids=uuids or [scan.uuid for scan in scans]))
    return received


def test_set_title_publishes_each_scans_own_metadata_value():
    scans = [make_titled_scan("s0", temp="1.5"), make_titled_scan("s1", temp="300")]

    received = set_title(scans, "temp")

    assert received[0].titles == [(UUID(value="s0"), "1.5"), (UUID(value="s1"), "300")]


def test_set_title_uses_the_text_as_typed_when_a_scan_has_no_such_field():
    received = set_title([make_titled_scan("s0", temp="1.5"), make_titled_scan("s1")], "temp")

    assert received[0].titles == [(UUID(value="s0"), "1.5"), (UUID(value="s1"), "temp")]


def test_set_title_leaves_out_a_scan_whose_field_is_blank():
    """SPICE files often carry an empty scan_title; the scan keeps its name rather than a blank label."""
    received = set_title([make_titled_scan("s0", scan_title="")], "scan_title")

    assert received[0].titles == []


def test_set_title_leaves_out_a_scan_no_longer_loaded():
    received = set_title([make_titled_scan("s0", temp="1.5")], "temp", uuids=[UUID(value="s0"), UUID(value="gone")])

    assert received[0].titles == [(UUID(value="s0"), "1.5")]


def test_a_blank_title_clears_every_title():
    received = set_title([make_titled_scan("s0", temp="1.5")], "  ")

    assert received[0].titles == []
