"""Tests for tavi.frontend.presenter.load_raw_scan_presenter."""

from unittest.mock import MagicMock, patch

from tavi.frontend.presenter.load_raw_scan_presenter import LoadRawScanPresenter
from tavi.library.data.scan import UUID
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.model_event import (
    AddFitEvent,
    AddPlotEvent,
    AddRawScanEvent,
    RemoveFitEvent,
    RemovePlotEvent,
    RemoveRawScanEvent,
)
from tavi.meta.event.type.presenter_event import FocusEvent


def _make_presenter():
    """Build a LoadRawScanPresenter with mock view and model."""
    view = MagicMock()
    model = MagicMock()
    with patch.object(LoadRawScanPresenter, 'init_view', lambda self: setattr(self, '_view', view)):
        presenter = LoadRawScanPresenter(model)
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


def test_init_inventory_is_empty():
    presenter, view, model = _make_presenter()

    assert presenter.inventory == {}


# ---------------------------------------------------------------------------
# update_treeview_data
# ---------------------------------------------------------------------------


def test_update_treeview_data_calls_add_raw_scan():
    presenter, view, model = _make_presenter()

    uuid = UUID(value="scan-001")
    event = AddRawScanEvent(uuid=uuid, friendly_name="My Scan", friendly_path="/exp1")
    presenter.update_treeview_data(event)

    view.add_raw_scan.assert_called_once_with(uuid, "My Scan", "/exp1")


def test_update_treeview_data_updates_inventory():
    presenter, view, model = _make_presenter()

    uuid = UUID(value="scan-002")
    event = AddRawScanEvent(uuid=uuid, friendly_name="Scan B", friendly_path="/exp2")
    presenter.update_treeview_data(event)

    assert presenter.inventory[uuid] == ("Scan B", "/exp2")


def test_update_treeview_data_accumulates_multiple_scans():
    presenter, view, model = _make_presenter()

    uuids = [UUID(value=f"u{i}") for i in range(3)]
    for i, uuid in enumerate(uuids):
        event = AddRawScanEvent(uuid=uuid, friendly_name=f"Scan{i}", friendly_path=f"/exp{i}")
        presenter.update_treeview_data(event)

    assert len(presenter.inventory) == 3
    for i, uuid in enumerate(uuids):
        assert uuid in presenter.inventory


def test_update_treeview_data_via_event_broker():
    presenter, view, model = _make_presenter()
    broker = EventBroker()

    uuid = UUID(value="broker-1")
    broker.publish(AddRawScanEvent(uuid=uuid, friendly_name="BrokerScan", friendly_path="/broker"))

    view.add_raw_scan.assert_called_once_with(uuid, "BrokerScan", "/broker")
    assert presenter.inventory[uuid] == ("BrokerScan", "/broker")


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


def test_update_plot_treeview_data_does_not_touch_raw_scan_inventory():
    presenter, view, model = _make_presenter()

    uuid = UUID(value="plot-002")
    presenter.update_plot_treeview_data(AddPlotEvent(uuid=uuid, friendly_name="run2_Plot", friendly_path=""))

    assert uuid not in presenter.inventory


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


def test_update_fit_treeview_data_does_not_touch_raw_scan_inventory():
    presenter, view, model = _make_presenter()

    uuid = UUID(value="fit-002")
    presenter.update_fit_treeview_data(AddFitEvent(uuid=uuid, friendly_name="run2_Fit", friendly_path=""))

    assert uuid not in presenter.inventory


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


def test_remove_treeview_data_drops_inventory_entry():
    presenter, view, model = _make_presenter()
    uuid = UUID(value="u1")
    presenter.update_treeview_data(AddRawScanEvent(uuid=uuid, friendly_name="scan", friendly_path="/exp"))

    presenter.remove_treeview_data(RemoveRawScanEvent(uuid=uuid))

    assert uuid not in presenter.inventory


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


def test_remove_treeview_data_tolerates_uuid_not_in_inventory():
    presenter, view, model = _make_presenter()

    presenter.remove_treeview_data(RemovePlotEvent(uuid=UUID(value="never-tracked")))

    view.remove_item.assert_called_once()
