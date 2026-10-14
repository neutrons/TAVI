"""Tests for FitHistory."""

from tavi.backend.model.fit_history import FitHistory
from tavi.library.data.fit_entry import FitMember, ParamField, PeakField
from tavi.library.data.plot import PlotSeries
from tavi.library.data.scan import UUID

FIT = UUID(value="fit-001")


def make_param(value=0) -> ParamField:
    return ParamField(value=str(value), fixed=False, minimum="", maximum="")


def make_member(range_max="10", uuid_val="scan-001") -> FitMember:
    return FitMember(
        series=PlotSeries(
            source_scan_uuid=UUID(value=uuid_val),
            scan_name="test_scan",
            normalized_by=None,
            x_name="qh",
            y_name="en",
            error_name="err",
        ),
        range_min="0",
        range_max=range_max,
        background="None",
        background_constant=make_param(0),
        peaks=[PeakField(shape="Gaussian", amplitude=make_param(1), center=make_param(0), fwhm=make_param(1))],
    )


def key(uuid_val="scan-001") -> tuple[UUID, UUID]:
    return (FIT, UUID(value=uuid_val))


def test_new_history_has_nothing_to_undo_or_redo():
    history = FitHistory()

    assert not history.can_undo(key())
    assert not history.can_redo(key())
    assert history.undo(FIT, make_member()) is None
    assert history.redo(FIT, make_member()) is None


def test_undo_returns_the_recorded_state_and_keeps_current_for_redo():
    history = FitHistory()
    history.record(FIT, make_member("1"))

    restored = history.undo(FIT, make_member("2"))

    assert restored.range_max == "1"
    assert not history.can_undo(key())
    assert history.redo(FIT, restored).range_max == "2"
    assert history.can_undo(key())


def test_undo_steps_back_in_reverse_save_order():
    history = FitHistory()
    history.record(FIT, make_member("1"))
    history.record(FIT, make_member("2"))

    assert history.undo(FIT, make_member("3")).range_max == "2"
    assert history.undo(FIT, make_member("2")).range_max == "1"


def test_recording_a_new_step_clears_redo():
    history = FitHistory()
    history.record(FIT, make_member("1"))
    history.undo(FIT, make_member("2"))

    history.record(FIT, make_member("1"))

    assert not history.can_redo(key())


def test_limit_drops_the_oldest_state():
    history = FitHistory(limit=2)
    for range_max in ("1", "2", "3"):
        history.record(FIT, make_member(range_max))

    assert history.undo(FIT, make_member("4")).range_max == "3"
    assert history.undo(FIT, make_member("3")).range_max == "2"
    assert history.undo(FIT, make_member("2")) is None


def test_members_of_one_fit_have_separate_histories():
    history = FitHistory()
    history.record(FIT, make_member(uuid_val="scan-001"))

    assert history.can_undo(key("scan-001"))
    assert not history.can_undo(key("scan-002"))


def test_forget_drops_by_fit_and_by_scan():
    history = FitHistory()
    history.record(FIT, make_member(uuid_val="scan-001"))
    history.record(FIT, make_member(uuid_val="scan-002"))
    other_fit = UUID(value="fit-002")
    history.record(other_fit, make_member(uuid_val="scan-003"))

    history.forget(scan_uuids={UUID(value="scan-001")})
    assert not history.can_undo(key("scan-001"))
    assert history.can_undo(key("scan-002"))

    history.forget(fit_uuids={other_fit})
    assert not history.can_undo((other_fit, UUID(value="scan-003")))
    assert history.can_undo(key("scan-002"))
