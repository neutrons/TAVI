"""Tests for tavi.frontend.presenter.load_raw_scan_presenter."""

from unittest.mock import MagicMock, call, patch

import pytest

from tavi.frontend.presenter.load_raw_scan_presenter import LoadRawScanPresenter
from tavi.library.data.scan import UUID
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.model_event import PlotAppendEvent, RawScanAppendEvent
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

    assert presenter.update_treeview_data in broker.registry[RawScanAppendEvent]


def test_init_registers_focus_event():
    presenter, view, model = _make_presenter()
    broker = EventBroker()

    assert presenter.print_selected in broker.registry[FocusEvent]


def test_init_registers_plot_append_event():
    presenter, view, model = _make_presenter()
    broker = EventBroker()

    assert presenter.update_plot_treeview_data in broker.registry[PlotAppendEvent]


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
    event = RawScanAppendEvent(uuid=uuid, friendly_name="My Scan", friendly_path="/exp1")
    presenter.update_treeview_data(event)

    view.add_raw_scan.assert_called_once_with(uuid, "My Scan", "/exp1")


def test_update_treeview_data_updates_inventory():
    presenter, view, model = _make_presenter()

    uuid = UUID(value="scan-002")
    event = RawScanAppendEvent(uuid=uuid, friendly_name="Scan B", friendly_path="/exp2")
    presenter.update_treeview_data(event)

    assert presenter.inventory[uuid] == ("Scan B", "/exp2")


def test_update_treeview_data_accumulates_multiple_scans():
    presenter, view, model = _make_presenter()

    uuids = [UUID(value=f"u{i}") for i in range(3)]
    for i, uuid in enumerate(uuids):
        event = RawScanAppendEvent(uuid=uuid, friendly_name=f"Scan{i}", friendly_path=f"/exp{i}")
        presenter.update_treeview_data(event)

    assert len(presenter.inventory) == 3
    for i, uuid in enumerate(uuids):
        assert uuid in presenter.inventory


def test_update_treeview_data_via_event_broker():
    presenter, view, model = _make_presenter()
    broker = EventBroker()

    uuid = UUID(value="broker-1")
    broker.publish(RawScanAppendEvent(uuid=uuid, friendly_name="BrokerScan", friendly_path="/broker"))

    view.add_raw_scan.assert_called_once_with(uuid, "BrokerScan", "/broker")
    assert presenter.inventory[uuid] == ("BrokerScan", "/broker")


# ---------------------------------------------------------------------------
# update_plot_treeview_data
# ---------------------------------------------------------------------------


def test_update_plot_treeview_data_calls_add_plot():
    presenter, view, model = _make_presenter()

    uuid = UUID(value="plot-001")
    event = PlotAppendEvent(uuid=uuid, friendly_name="run1_Plot", friendly_path="")
    presenter.update_plot_treeview_data(event)

    view.add_plot.assert_called_once_with(uuid, "run1_Plot", "")


def test_update_plot_treeview_data_via_event_broker():
    presenter, view, model = _make_presenter()
    broker = EventBroker()

    uuid = UUID(value="plot-broker-1")
    broker.publish(PlotAppendEvent(uuid=uuid, friendly_name="BrokerPlot", friendly_path=""))

    view.add_plot.assert_called_once_with(uuid, "BrokerPlot", "")


def test_update_plot_treeview_data_does_not_touch_raw_scan_inventory():
    presenter, view, model = _make_presenter()

    uuid = UUID(value="plot-002")
    presenter.update_plot_treeview_data(PlotAppendEvent(uuid=uuid, friendly_name="run2_Plot", friendly_path=""))

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
# print_selected
# ---------------------------------------------------------------------------


def test_print_selected_does_not_raise(capsys):
    presenter, view, model = _make_presenter()

    event = FocusEvent(ids=[UUID(value="p1"), UUID(value="p2")])
    presenter.print_selected(event)

    captured = capsys.readouterr()
    assert "p1" in captured.out
    assert "p2" in captured.out
