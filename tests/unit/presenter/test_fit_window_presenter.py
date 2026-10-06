"""Tests for FitWindowPresenter."""

import pytest

from tavi.frontend.presenter.fit_window_presenter import FitWindowPresenter
from tavi.frontend.view.fit_window_view import FitWindow, FitWindowsView
from tavi.library.data.fit_entry import (
    FitCurve,
    FitData,
    FitMember,
    FitOutcome,
    FitResultSummary,
    ParamField,
    PeakField,
    PeakResult,
)
from tavi.library.data.plot import PlotSeries
from tavi.library.data.scan import UUID
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.model_event import FitRemoveEvent, RawScanRemoveEvent
from tavi.meta.event.type.presenter_event import SaveFitEvent, SyncFitEvent


def make_series(uuid_val="scan-001") -> PlotSeries:
    return PlotSeries(
        source_scan_uuid=UUID(value=uuid_val),
        scan_name="sample alignment",
        friendly_name=uuid_val,
        normalized_by=None,
        x_name="qh",
        y_name="detector",
        error_name="err",
    )


def make_param(value=0) -> ParamField:
    return ParamField(value=str(value), fixed=False, minimum="", maximum="")


def make_result(reduced_chi_squared=0.5) -> FitResultSummary:
    return FitResultSummary(
        reduced_chi_squared=reduced_chi_squared,
        peaks=[PeakResult(amplitude=1.0, amplitude_err=0.1, center=0.0, center_err=None, fwhm=1.0, fwhm_err=0.2)],
    )


def make_member(uuid_val="scan-001", reduced_chi_squared=0.5) -> FitMember:
    return FitMember(
        series=make_series(uuid_val),
        range_min="0",
        range_max="10",
        background="Linear",
        background_constant=make_param(0),
        peaks=[PeakField(shape="Gaussian", amplitude=make_param(1), center=make_param(0), fwhm=make_param(1))],
        result=make_result(reduced_chi_squared),
    )


def make_outcome(uuid_val="scan-001", reduced_chi_squared=0.5) -> FitOutcome:
    return FitOutcome(
        member=make_member(uuid_val, reduced_chi_squared),
        curve=FitCurve(source_scan_uuid=UUID(value=uuid_val), scan_name=uuid_val, x=[1.0, 2.0], best_fit=[1.1, 1.9]),
        data=FitData(x=[1.0, 2.0], y=[1.0, 2.0], err=[0.1, 0.1]),
    )


def make_sync_fit_event(uuid_vals=("scan-001", "scan-002"), fit_uuid="fit-001", reduced_chi_squared=0.5):
    return SyncFitEvent(
        fit_uuid=UUID(value=fit_uuid),
        outcomes=[make_outcome(uuid_val, reduced_chi_squared) for uuid_val in uuid_vals],
    )


def publish_performed_fit(sync: SyncFitEvent) -> None:
    """Publish what FitModel.perform_fit does: save the members, then sync the UI to them."""
    members = [outcome.member for outcome in sync.outcomes]
    EventBroker().publish(SaveFitEvent(fit_uuid=sync.fit_uuid, members=members))
    EventBroker().publish(sync)


def perform_fit(**kwargs) -> None:
    publish_performed_fit(make_sync_fit_event(**kwargs))


def recompute_fit(**kwargs) -> None:
    """Publish what a tree recompute does: sync the UI only, never save."""
    EventBroker().publish(make_sync_fit_event(**kwargs))


def key(fit_uuid="fit-001", uuid_val="scan-001") -> tuple[UUID, UUID]:
    return (UUID(value=fit_uuid), UUID(value=uuid_val))


@pytest.fixture
def presenter(qtbot):
    p = FitWindowPresenter()
    yield p
    for window in list(p._view.windows.values()):
        window.close()


# ---------------------------------------------------------------------------
# __init__ - wiring
# ---------------------------------------------------------------------------


def test_init_view_is_fit_windows_view(presenter):
    assert isinstance(presenter._view, FitWindowsView)


def test_init_starts_with_no_windows_open(presenter):
    assert presenter.open_windows() == set()


def test_init_registers_save_sync_and_removal_events(presenter):
    broker = EventBroker()
    assert presenter.handle_save_fit in broker.registry[SaveFitEvent]
    assert presenter.handle_sync_fit in broker.registry[SyncFitEvent]
    assert presenter.handle_fit_removed in broker.registry[FitRemoveEvent]
    assert presenter.handle_raw_scan_removed in broker.registry[RawScanRemoveEvent]


# ---------------------------------------------------------------------------
# handle_save_fit / handle_sync_fit - opening and refreshing
# ---------------------------------------------------------------------------


def test_multi_outcome_fit_opens_one_window_per_member(presenter):
    perform_fit()

    assert presenter.open_windows() == {key(uuid_val="scan-001"), key(uuid_val="scan-002")}
    assert set(presenter._view.windows) == presenter.open_windows()
    assert all(isinstance(window, FitWindow) for window in presenter._view.windows.values())


def test_opened_window_shows_its_members_result(presenter):
    perform_fit(reduced_chi_squared=0.25)

    window = presenter._view.windows[key(uuid_val="scan-002")]
    assert window.windowTitle() == "Fit - scan-002"
    assert window.status_label.text() == "χ2 = 0.25"
    assert window.param_table.rowCount() == 3


def test_single_outcome_fit_does_not_open_a_window(presenter):
    """A fit of one scan is already on the main plotter."""
    perform_fit(uuid_vals=("scan-001",))

    assert presenter.open_windows() == set()
    assert presenter._view.windows == {}


def test_single_outcome_refit_refreshes_an_open_window(presenter):
    perform_fit()
    window = presenter._view.windows[key(uuid_val="scan-001")]

    perform_fit(uuid_vals=("scan-001",), reduced_chi_squared=0.75)

    assert presenter._view.windows[key(uuid_val="scan-001")] is window
    assert window.status_label.text() == "χ2 = 0.75"


def test_rerunning_a_multi_outcome_fit_reuses_its_windows(presenter):
    perform_fit()
    windows = dict(presenter._view.windows)

    perform_fit(reduced_chi_squared=0.75)

    assert presenter._view.windows == windows
    assert windows[key(uuid_val="scan-002")].status_label.text() == "χ2 = 0.75"


def test_failed_member_window_says_the_fit_did_not_run(presenter):
    failed = FitOutcome(member=make_member("scan-002").model_copy(update={"result": None}))
    publish_performed_fit(SyncFitEvent(fit_uuid=UUID(value="fit-001"), outcomes=[make_outcome(), failed]))

    window = presenter._view.windows[key(uuid_val="scan-002")]
    assert "did not run" in window.status_label.text()
    assert window.param_table.rowCount() == 0


# ---------------------------------------------------------------------------
# user closing a window
# ---------------------------------------------------------------------------


def test_user_closing_a_window_removes_it_from_the_registry(presenter):
    perform_fit()

    presenter._view.windows[key(uuid_val="scan-001")].close()

    assert presenter.open_windows() == {key(uuid_val="scan-002")}
    assert set(presenter._view.windows) == {key(uuid_val="scan-002")}


def test_single_outcome_refit_does_not_reopen_a_closed_window(presenter):
    perform_fit()
    presenter._view.windows[key(uuid_val="scan-001")].close()

    perform_fit(uuid_vals=("scan-001",))

    assert key(uuid_val="scan-001") not in presenter.open_windows()
    assert key(uuid_val="scan-001") not in presenter._view.windows


def test_reselecting_a_fit_from_the_tree_opens_no_windows(presenter):
    """A recompute (the tree path) shows the fit on the main plotter only."""
    recompute_fit()

    assert presenter.open_windows() == set()
    assert presenter._view.windows == {}


def test_recompute_leaves_open_windows_untouched(presenter):
    """Browsing a fit from the tree only shows it on the main plotter, even if its windows are open."""
    perform_fit()
    window = presenter._view.windows[key(uuid_val="scan-001")]

    recompute_fit(reduced_chi_squared=0.75)

    assert window.status_label.text() == "χ2 = 0.5"
    assert presenter.open_windows() == {key(uuid_val="scan-001"), key(uuid_val="scan-002")}


def test_close_all_button_closes_every_fit_window(presenter):
    perform_fit()
    perform_fit(fit_uuid="fit-002")
    assert len(presenter.open_windows()) == 4

    presenter._view.windows[key(uuid_val="scan-001")].close_all_btn.click()

    assert presenter.open_windows() == set()
    assert presenter._view.windows == {}


def test_multi_outcome_refit_reopens_a_closed_window(presenter):
    """Rerunning the whole sequential fit shows every member again."""
    perform_fit()
    presenter._view.windows[key(uuid_val="scan-001")].close()

    perform_fit()

    assert key(uuid_val="scan-001") in presenter.open_windows()
    assert key(uuid_val="scan-001") in presenter._view.windows


# ---------------------------------------------------------------------------
# removals
# ---------------------------------------------------------------------------


def test_fit_removed_closes_its_windows_only(presenter):
    perform_fit(fit_uuid="fit-001")
    perform_fit(fit_uuid="fit-002")

    EventBroker().publish(FitRemoveEvent(uuid=UUID(value="fit-001")))

    assert presenter.open_windows() == {key("fit-002", "scan-001"), key("fit-002", "scan-002")}
    assert set(presenter._view.windows) == presenter.open_windows()


def test_raw_scan_removed_closes_every_window_showing_that_scan(presenter):
    perform_fit(fit_uuid="fit-001")
    perform_fit(fit_uuid="fit-002")

    EventBroker().publish(RawScanRemoveEvent(uuid=UUID(value="scan-001")))

    assert presenter.open_windows() == {key("fit-001", "scan-002"), key("fit-002", "scan-002")}
    assert set(presenter._view.windows) == presenter.open_windows()


def test_removal_matching_no_window_is_a_noop(presenter):
    perform_fit()

    EventBroker().publish(FitRemoveEvent(uuid=UUID(value="fit-999")))
    EventBroker().publish(RawScanRemoveEvent(uuid=UUID(value="scan-999")))

    assert len(presenter.open_windows()) == 2


def test_open_windows_returns_a_copy(presenter):
    """Callers must not be able to mutate the registry through the returned set."""
    perform_fit()

    presenter.open_windows().clear()

    assert len(presenter.open_windows()) == 2
