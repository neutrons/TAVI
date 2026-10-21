"""Tests for tavi.frontend.presenter.project_presenter."""

from unittest.mock import MagicMock, patch

from tavi.frontend.presenter.project_presenter import ProjectPresenter
from tavi.library.data.scan import UUID
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.model_event import (
    AddFitEvent,
    AddPlotEvent,
    AddRawScanEvent,
    RemoveFitEvent,
    RemovePlotEvent,
    ClearRawScansEvent,
    RemoveRawScanEvent,
)
from tavi.meta.event.type.presenter_event import FocusEvent


def _make_presenter():
    """Build a ProjectPresenter with mock view and model."""
    view = MagicMock()
    model = MagicMock()
    with patch.object(ProjectPresenter, 'init_view', lambda self: setattr(self, '_view', view)):
        presenter = ProjectPresenter(model)
    return presenter, view, model


# ---------------------------------------------------------------------------
# __init__ — wiring
# ---------------------------------------------------------------------------


def test_init_registers_raw_scan_append_event():
    presenter, view, model = _make_presenter()
    broker = EventBroker()

    assert presenter.update_treeview_data in broker.registry[AddRawScanEvent]


def test_init_registers_plot_append_event():
    presenter, view, model = _make_presenter()
    broker = EventBroker()

    assert presenter.update_plot_treeview_data in broker.registry[AddPlotEvent]


def test_init_registers_fit_append_event():
    presenter, view, model = _make_presenter()
    broker = EventBroker()

    assert presenter.update_fit_treeview_data in broker.registry[AddFitEvent]


def test_init_hooks_up_select_signal():
    presenter, view, model = _make_presenter()

    view.hookup_select_signal.assert_called_once_with(presenter.handle_selection_event)


# ---------------------------------------------------------------------------
# update_treeview_data
# ---------------------------------------------------------------------------


def test_update_treeview_data_calls_add_raw_scan():
    presenter, view, model = _make_presenter()

    uuid = UUID(value="scan-001")
    event = AddRawScanEvent(uuid=uuid, friendly_name="My Scan", friendly_path="/exp1")
    presenter.update_treeview_data(event)

    view.add_raw_scan.assert_called_once_with(uuid, "My Scan", "/exp1")


def test_update_treeview_data_via_event_broker():
    presenter, view, model = _make_presenter()
    broker = EventBroker()

    uuid = UUID(value="broker-1")
    broker.publish(AddRawScanEvent(uuid=uuid, friendly_name="BrokerScan", friendly_path="/broker"))

    view.add_raw_scan.assert_called_once_with(uuid, "BrokerScan", "/broker")


# ---------------------------------------------------------------------------
# update_plot_treeview_data
# ---------------------------------------------------------------------------


def test_update_plot_treeview_data_calls_add_plot():
    presenter, view, model = _make_presenter()

    uuid = UUID(value="plot-001")
    event = AddPlotEvent(uuid=uuid, friendly_name="run1_Plot", friendly_path="")
    presenter.update_plot_treeview_data(event)

    view.add_plot.assert_called_once_with(uuid, "run1_Plot", "")


def test_update_plot_treeview_data_via_event_broker():
    presenter, view, model = _make_presenter()
    broker = EventBroker()

    uuid = UUID(value="plot-broker-1")
    broker.publish(AddPlotEvent(uuid=uuid, friendly_name="BrokerPlot", friendly_path=""))

    view.add_plot.assert_called_once_with(uuid, "BrokerPlot", "")


# ---------------------------------------------------------------------------
# update_fit_treeview_data
# ---------------------------------------------------------------------------


def test_update_fit_treeview_data_calls_add_fit():
    presenter, view, model = _make_presenter()

    uuid = UUID(value="fit-001")
    event = AddFitEvent(uuid=uuid, friendly_name="run1_Fit", friendly_path="")
    presenter.update_fit_treeview_data(event)

    view.add_fit.assert_called_once_with(uuid, "run1_Fit", "")


def test_update_fit_treeview_data_via_event_broker():
    presenter, view, model = _make_presenter()
    broker = EventBroker()

    uuid = UUID(value="fit-broker-1")
    broker.publish(AddFitEvent(uuid=uuid, friendly_name="BrokerFit", friendly_path=""))

    view.add_fit.assert_called_once_with(uuid, "BrokerFit", "")


# ---------------------------------------------------------------------------
# handle_selection_event
# ---------------------------------------------------------------------------


def test_handle_selection_event_publishes_focus_event():
    presenter, view, model = _make_presenter()
    broker = EventBroker()

    uuid = UUID(value="focus-1")
    view.get_selected_items.return_value = [uuid]

    received: list[FocusEvent] = []
    broker.register(FocusEvent, received.append)

    presenter.handle_selection_event()

    assert len(received) == 1
    assert received[0].ids == [uuid]


def test_handle_selection_event_uses_view_selection():
    presenter, view, model = _make_presenter()

    uuids = [UUID(value="a"), UUID(value="b")]
    view.get_selected_items.return_value = uuids

    received: list[FocusEvent] = []
    EventBroker().register(FocusEvent, received.append)

    presenter.handle_selection_event()

    view.get_selected_items.assert_called_once()
    assert received[0].ids == uuids


def test_handle_selection_event_empty_selection():
    presenter, view, model = _make_presenter()

    view.get_selected_items.return_value = []

    received: list[FocusEvent] = []
    EventBroker().register(FocusEvent, received.append)

    presenter.handle_selection_event()

    assert received[0].ids == []


# ---------------------------------------------------------------------------
# removal
# ---------------------------------------------------------------------------


def test_init_registers_raw_scan_remove_event():
    presenter, view, model = _make_presenter()
    broker = EventBroker()

    assert presenter.remove_treeview_data in broker.registry[RemoveRawScanEvent]


def test_init_registers_plot_remove_event():
    presenter, view, model = _make_presenter()
    broker = EventBroker()

    assert presenter.remove_treeview_data in broker.registry[RemovePlotEvent]


def test_init_registers_fit_remove_event():
    presenter, view, model = _make_presenter()
    broker = EventBroker()

    assert presenter.remove_treeview_data in broker.registry[RemoveFitEvent]


def test_init_hooks_up_remove_signal():
    presenter, view, model = _make_presenter()

    view.hookup_remove_signal.assert_called_once_with(presenter.handle_remove_request)


def test_handle_remove_request_forwards_to_model():
    presenter, view, model = _make_presenter()
    uuids = [UUID(value="u1"), UUID(value="u2")]

    presenter.handle_remove_request(uuids)

    model.remove_items.assert_called_once_with(uuids)


def test_handle_remove_request_does_not_touch_view():
    """The tree must wait for the model's remove event, not delete rows on request."""
    presenter, view, model = _make_presenter()

    presenter.handle_remove_request([UUID(value="u1")])

    view.remove_item.assert_not_called()


def test_remove_treeview_data_removes_from_view():
    presenter, view, model = _make_presenter()
    uuid = UUID(value="u1")

    presenter.remove_treeview_data(RemoveRawScanEvent(uuid=uuid))

    view.remove_item.assert_called_once_with(uuid)


def test_remove_treeview_data_handles_plot_remove_event():
    presenter, view, model = _make_presenter()
    uuid = UUID(value="p1")

    presenter.remove_treeview_data(RemovePlotEvent(uuid=uuid))

    view.remove_item.assert_called_once_with(uuid)


def test_remove_treeview_data_handles_fit_remove_event():
    presenter, view, model = _make_presenter()
    uuid = UUID(value="f1")

    presenter.remove_treeview_data(RemoveFitEvent(uuid=uuid))

    view.remove_item.assert_called_once_with(uuid)


# ---------------------------------------------------------------------------
# clear_raw_scan_treeview_data
# ---------------------------------------------------------------------------


def test_init_registers_clear_raw_scans_event():
    presenter, view, model = _make_presenter()

    assert presenter.clear_raw_scan_treeview_data in EventBroker().registry[ClearRawScansEvent]


def test_clear_raw_scans_event_clears_only_raw_scans_in_view():
    presenter, view, model = _make_presenter()

    EventBroker().publish(ClearRawScansEvent())

    view.clear_raw_scans.assert_called_once_with()
    view.remove_item.assert_not_called()
