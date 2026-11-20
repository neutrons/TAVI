"""Tests for tavi.frontend.presenter.project_presenter."""

from unittest.mock import MagicMock, patch

from qtpy.QtCore import QItemSelectionModel, QPoint
from qtpy.QtWidgets import QMenu

from tavi.backend.model.tavi_project_model import TaviProjectModel
from tavi.frontend.presenter.project_presenter import ProjectPresenter
from tavi.library.data.fit_entry import FitEntry, FitMember, ParamField, PeakField
from tavi.library.data.plot import Plot, PlotSeries
from tavi.library.data.scan import UUID, ProcessedScan, Provenance, RawScan, ScanData, ScanMetadata, TaviMetadata
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.model_event import (
    AddFitEvent,
    AddPlotEvent,
    AddProcessedScanEvent,
    AddRawScanEvent,
    RemoveFitEvent,
    RemovePlotEvent,
    RemoveProcessedScanEvent,
    ClearRawScansEvent,
    RemoveRawScanEvent,
    SyncItemNameEvent,
)
from tavi.meta.event.type.presenter_event import ClearFocusEvent, FocusEvent, GroupScansEvent


def make_raw_scan(uuid_val):
    return RawScan(
        uuid=UUID(value=uuid_val),
        data=ScanData(data={"qh": [1.0, 2.0], "en": [3.0, 4.0]}),
        metadata=ScanMetadata(),
        tavimeta=TaviMetadata(default_axis=("qh", "en"), friendly_name=uuid_val, friendly_path="/exp1"),
        prov=Provenance(raw_file="scan.dat", contributing_scans={UUID(value=uuid_val): 1}),
    )


def make_processed_scan(uuid_val):
    return ProcessedScan(
        uuid=UUID(value=uuid_val),
        data=ScanData(data={"qh": [1.0, 2.0], "en": [3.0, 4.0]}),
        metadata=ScanMetadata(),
        tavimeta=TaviMetadata(default_axis=("qh", "en"), friendly_name=uuid_val, friendly_path=""),
        prov=Provenance(raw_file="", contributing_scans={UUID(value="scan-001"): 1, UUID(value="scan-002"): 1}),
    )


def make_series(source_uuid):
    return PlotSeries(
        source_scan_uuid=UUID(value=source_uuid),
        scan_name=source_uuid,
        normalized_by=None,
        x_name="qh",
        y_name="en",
        error_name="err",
    )


def make_plot(uuid_val, source_uuid):
    return Plot(uuid=UUID(value=uuid_val), series=[make_series(source_uuid)])


def make_fit_entry(uuid_val, source_uuid):
    param = ParamField(value="0", fixed=False, minimum="", maximum="")
    member = FitMember(
        series=make_series(source_uuid),
        range_min="0",
        range_max="10",
        background="None",
        background_constant=param,
        peaks=[PeakField(shape="Gaussian", amplitude=param, center=param, fwhm=param)],
    )
    return FitEntry(uuid=UUID(value=uuid_val), members=[member])


def make_project_model():
    with (
        patch("tavi.backend.model.tavi_project_model.Resource"),
        patch("tavi.backend.model.tavi_project_model.RawScanLoadController"),
    ):
        return TaviProjectModel(MagicMock())


def _make_presenter():
    """Build a ProjectPresenter with mock view and model."""
    view = MagicMock()
    model = MagicMock()
    with patch.object(ProjectPresenter, "init_view", lambda self: setattr(self, "_view", view)):
        presenter = ProjectPresenter(model)
    return presenter, view, model


def make_scan(uuid_val="scan-001", name="scan0001", columns=None) -> RawScan:
    if columns is None:
        columns = {"qh": [1.0, 2.0], "detector": [3.0, 4.0], "mcu": [1.0, 1.0]}
    return RawScan(
        uuid=UUID(value=uuid_val),
        data=ScanData(data=columns),
        metadata=ScanMetadata(),
        tavimeta=TaviMetadata(default_axis=("qh", "detector"), friendly_name=name, friendly_path="/exp"),
        prov=Provenance(raw_file=f"{name}.dat", contributing_scans={UUID(value=uuid_val): 1}),
    )



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


def test_init_registers_processed_scan_remove_event():
    presenter, view, model = _make_presenter()
    broker = EventBroker()

    assert presenter.remove_treeview_data in broker.registry[RemoveProcessedScanEvent]


def test_processed_scan_remove_event_drops_the_tree_row():
    presenter, view, model = _make_presenter()
    uuid = UUID(value="hist-001")

    EventBroker().publish(RemoveProcessedScanEvent(uuid=uuid))

    view.remove_item.assert_called_once_with(uuid)


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
    scan = make_scan(uuid_val="u1")

    presenter.remove_treeview_data(RemoveRawScanEvent(uuid=scan.uuid, raw_scan=scan))

    view.remove_item.assert_called_once_with(scan.uuid)


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


# ---------------------------------------------------------------------------
# renaming
# ---------------------------------------------------------------------------


def test_init_registers_sync_item_name_event():
    presenter, view, model = _make_presenter()

    assert presenter.rename_treeview_data in EventBroker().registry[SyncItemNameEvent]


def test_init_hooks_up_rename_signal():
    presenter, view, model = _make_presenter()

    view.hookup_rename_signal.assert_called_once_with(presenter.handle_rename_request)


def test_handle_rename_request_forwards_to_model():
    presenter, view, model = _make_presenter()
    uuid = UUID(value="p1")

    presenter.handle_rename_request(uuid, "Bragg peak")

    model.rename_item.assert_called_once_with(uuid, "Bragg peak")


def test_handle_rename_request_does_not_touch_view():
    """The tree must wait for the model's sync event, not relabel the row on request."""
    presenter, view, model = _make_presenter()

    presenter.handle_rename_request(UUID(value="p1"), "Bragg peak")

    view.rename_item.assert_not_called()


def test_rename_treeview_data_relabels_in_view():
    presenter, view, model = _make_presenter()
    uuid = UUID(value="p1")

    presenter.rename_treeview_data(SyncItemNameEvent(uuid=uuid, friendly_name="Bragg peak"))

    view.rename_item.assert_called_once_with(uuid, "Bragg peak")


def test_rename_treeview_data_via_event_broker():
    presenter, view, model = _make_presenter()
    uuid = UUID(value="f1")

    EventBroker().publish(SyncItemNameEvent(uuid=uuid, friendly_name="Gaussian over run1"))

    view.rename_item.assert_called_once_with(uuid, "Gaussian over run1")


def test_add_processed_scan_event_lists_the_combined_scan():
    """A combination the model announces is listed in the tree, under its Combined root."""
    _presenter, view, _model = _make_presenter()
    uuid = UUID(value="combined-1")

    EventBroker().publish(AddProcessedScanEvent(uuid=uuid, friendly_name="scan0001+scan0002", friendly_path=""))

    view.add_combined.assert_called_once_with(uuid, "scan0001+scan0002", "")


# ---------------------------------------------------------------------------
# Add to Group Fit
# ---------------------------------------------------------------------------


def test_init_hooks_up_group_fit_signal():
    presenter, view, model = _make_presenter()

    view.hookup_group_fit_signal.assert_called_once_with(presenter.handle_group_fit_request)


def test_handle_group_fit_request_publishes_the_scans_in_order():
    presenter, view, model = _make_presenter()
    received = []
    EventBroker().register(GroupScansEvent, received.append)
    scans = [(UUID(value="s2"), "scan2"), (UUID(value="s1"), "scan1")]

    presenter.handle_group_fit_request(scans)

    assert [event.scans for event in received] == [scans]


# ---------------------------------------------------------------------------
# Removing a multi-selection, end to end
# ---------------------------------------------------------------------------


def test_removing_a_mixed_multi_selection_does_not_refocus_purged_items(qtbot, monkeypatch):
    """
    The model purges a whole batch before the tree drops its first row. Were the rows still
    selected, each removal would re-focus batch-mates TaviData no longer holds, and crash.
    """
    model = make_project_model()
    presenter = ProjectPresenter(model)
    tree = presenter._view.tree_widget
    qtbot.addWidget(presenter._view)

    raw = [make_raw_scan("scan-001"), make_raw_scan("scan-002")]
    combined = make_processed_scan("combined-001")
    plot = make_plot("plot-001", "scan-001")
    fit = make_fit_entry("fit-001", "scan-002")
    for scan in raw:
        model.tavi_data.raw_scans[scan.uuid] = scan
        tree.add_raw_scan(scan.uuid, scan.uuid.value, "/exp1")
    model.tavi_data.processed_scans[combined.uuid] = combined
    tree.add_combined(combined.uuid, "combined", "")
    model.tavi_data.plots[plot.uuid] = plot
    tree.add_plot(plot.uuid, "plot", "")
    model.tavi_data.fits[fit.uuid] = fit
    tree.add_fit(fit.uuid, "fit", "")

    selected = [scan.uuid for scan in raw] + [combined.uuid, plot.uuid, fit.uuid]
    for uuid in selected:
        tree.treeView.selectionModel().select(tree.uuid_map[uuid].index(), QItemSelectionModel.Select)

    cleared = []
    EventBroker().register(ClearFocusEvent, cleared.append)
    monkeypatch.setattr(QMenu, "exec", lambda self, *a, **k: next(a for a in self.actions() if a.text() == "Remove Item"))
    monkeypatch.setattr(tree.treeView, "indexAt", lambda _pos: tree.uuid_map[selected[0]].index())
    tree.show_context_menu(QPoint(0, 0))

    assert tree.uuid_map == {}
    assert model.tavi_data.raw_scans == {}
    assert model.tavi_data.processed_scans == {}
    assert model.tavi_data.plots == {}
    assert model.tavi_data.fits == {}
    assert len(cleared) == 1
