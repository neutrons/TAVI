"""Tests for TaviProjectModel."""

from unittest.mock import MagicMock, patch

import pytest

from tavi.backend.model.tavi_project_model import TaviProjectModel
from tavi.library.data.model_response import ModelResponse, ResponseCode
from tavi.library.data.plot import Plot, PlotSeries
from tavi.library.data.scan import UUID, Provenance, RawScan, ScanData, ScanMetadata, TaviMetadata
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.model_event import PlotAppendEvent, RawScanAppendEvent, SyncRecentProjects
from tavi.meta.event.type.presenter_event import (
    ActivePlotChangedEvent,
    DownstreamReadyEvent,
    FocusActivePlotEvent,
    FocusEvent,
    PlotFocusEvent,
    RawScanFocusEvent,
    SavePlotEvent,
)

SETTINGS_YAML = "TAVI:\n  recent:\n    projects:\n      - /path/to/project1\n      - /path/to/project2\n"


def make_filestore():
    fs = MagicMock()
    fs.read_user_data_file.return_value = SETTINGS_YAML
    return fs


def make_raw_scan(uuid_val="scan-001"):
    return RawScan(
        uuid=UUID(value=uuid_val),
        data=ScanData(data={"qh": [1.0, 2.0], "en": [3.0, 4.0]}),
        metadata=ScanMetadata(),
        tavimeta=TaviMetadata(
            default_axis=("qh", "en"),
            normalization=("monitor", 1.0),
            friendly_name="test_scan",
            friendly_path="/exp1",
        ),
        prov=Provenance(raw_file="scan.dat", contributing_scans={UUID(value=uuid_val): 1}),
    )


def make_plot(uuid_val="plot-001", series=None):
    if series is None:
        series = [
            PlotSeries(
                source_scan_uuid=UUID(value="scan-001"),
                scan_name="test_plot",
                normalized_by="monitor",
                x_name="qh",
                y_name="en",
                error_name="err",
            )
        ]
    return Plot(uuid=UUID(value=uuid_val), series=series)


def make_series(scan_name, uuid_val="scan-001"):
    return PlotSeries(
        source_scan_uuid=UUID(value=uuid_val),
        scan_name=scan_name,
        normalized_by=None,
        x_name="qh",
        y_name="en",
        error_name="err",
    )


@pytest.fixture(autouse=True)
def patch_resource():
    with patch("tavi.backend.model.tavi_project_model.Resource"):
        yield


@pytest.fixture
def model():
    with patch("tavi.backend.model.tavi_project_model.RawScanLoadController"):
        m = TaviProjectModel(make_filestore())
    return m


# ---------------------------------------------------------------------------
# __init__ — event wiring
# ---------------------------------------------------------------------------


def test_init_registers_downstream_ready_handler(model):
    broker = EventBroker()
    assert model.sync_on_ready in broker.registry[DownstreamReadyEvent]


def test_init_registers_focus_event_handler(model):
    broker = EventBroker()
    assert model._handle_focus_event in broker.registry[FocusEvent]


def test_init_tavi_data_empty(model):
    assert model.tavi_data.raw_scans == {}
    assert model.tavi_data.plots == {}


# ---------------------------------------------------------------------------
# get_plots_handle
# ---------------------------------------------------------------------------


def test_get_plots_handle_returns_same_dict(model):
    handle = model.get_plots_handle()
    assert handle is model.tavi_data.plots


def test_get_plots_handle_reflects_mutations(model):
    plot = make_plot()
    handle = model.get_plots_handle()
    model.tavi_data.plots[plot.uuid] = plot
    assert plot.uuid in handle


# ---------------------------------------------------------------------------
# load_raw_scan_from_folder
# ---------------------------------------------------------------------------


def test_load_raw_scan_from_folder_returns_ok(model):
    model.raw_scan_load_controller = MagicMock()
    model.raw_scan_load_controller.load_folder.return_value = []
    result = model.load_raw_scan_from_folder("/some/folder")
    assert result.code == ResponseCode.OK


def test_load_raw_scan_from_folder_stores_scan(model):
    scan = make_raw_scan()
    model.raw_scan_load_controller = MagicMock()
    model.raw_scan_load_controller.load_folder.return_value = [scan]

    model.load_raw_scan_from_folder("/some/folder")

    assert scan.uuid in model.tavi_data.raw_scans


def test_load_raw_scan_from_folder_publishes_append_event(model):
    scan = make_raw_scan()
    model.raw_scan_load_controller = MagicMock()
    model.raw_scan_load_controller.load_folder.return_value = [scan]

    received = []
    EventBroker().register(RawScanAppendEvent, received.append)
    model.load_raw_scan_from_folder("/some/folder")

    assert len(received) == 1
    assert received[0].uuid == scan.uuid
    assert received[0].friendly_name == scan.tavimeta.friendly_name


def test_load_raw_scan_from_folder_multiple_scans_all_stored(model):
    scans = [make_raw_scan(f"scan-{i:03d}") for i in range(3)]
    model.raw_scan_load_controller = MagicMock()
    model.raw_scan_load_controller.load_folder.return_value = scans

    model.load_raw_scan_from_folder("/some/folder")

    assert len(model.tavi_data.raw_scans) == 3
    for scan in scans:
        assert scan.uuid in model.tavi_data.raw_scans


def test_load_raw_scan_from_folder_multiple_scans_each_publishes_event(model):
    scans = [make_raw_scan(f"scan-{i:03d}") for i in range(3)]
    model.raw_scan_load_controller = MagicMock()
    model.raw_scan_load_controller.load_folder.return_value = scans

    received = []
    EventBroker().register(RawScanAppendEvent, received.append)
    model.load_raw_scan_from_folder("/some/folder")

    assert len(received) == 3


def test_load_raw_scan_from_folder_empty_folder(model):
    model.raw_scan_load_controller = MagicMock()
    model.raw_scan_load_controller.load_folder.return_value = []

    received = []
    EventBroker().register(RawScanAppendEvent, received.append)
    model.load_raw_scan_from_folder("/empty/folder")

    assert len(received) == 0
    assert model.tavi_data.raw_scans == {}


# ---------------------------------------------------------------------------
# sync_on_ready
# ---------------------------------------------------------------------------


def test_sync_on_ready_calls_emit_sync_recent_projects(model):
    model.emit_sync_recent_projects = MagicMock()
    model.sync_on_ready(DownstreamReadyEvent())
    model.emit_sync_recent_projects.assert_called_once()


# ---------------------------------------------------------------------------
# emit_sync_recent_projects
# ---------------------------------------------------------------------------


def test_emit_sync_recent_projects_publishes_event(model):
    received = []
    EventBroker().register(SyncRecentProjects, received.append)

    model.emit_sync_recent_projects()

    assert len(received) == 1
    assert isinstance(received[0], SyncRecentProjects)


def test_emit_sync_recent_projects_includes_projects_from_settings(model):
    received = []
    EventBroker().register(SyncRecentProjects, received.append)

    model.emit_sync_recent_projects()

    assert "/path/to/project1" in received[0].recent_projects
    assert "/path/to/project2" in received[0].recent_projects


def test_emit_sync_recent_projects_reads_from_filestore(model):
    model.emit_sync_recent_projects()
    model.filestore.read_user_data_file.assert_called_with("settings.yaml")


def test_downstream_ready_event_triggers_sync(model):
    received = []
    EventBroker().register(SyncRecentProjects, received.append)

    EventBroker().publish(DownstreamReadyEvent())

    assert len(received) == 1


# ---------------------------------------------------------------------------
# _handle_focus_event — raw scan routing
# ---------------------------------------------------------------------------


def test_handle_focus_event_raw_scan_publishes_raw_scan_focus_event(model):
    scan = make_raw_scan()
    model.tavi_data.raw_scans[scan.uuid] = scan

    received = []
    EventBroker().register(RawScanFocusEvent, received.append)
    EventBroker().publish(FocusEvent(ids=[scan.uuid]))

    assert len(received) == 1
    assert received[0].scans[0].uuid == scan.uuid


def test_handle_focus_event_raw_scan_does_not_publish_plot_event(model):
    scan = make_raw_scan()
    model.tavi_data.raw_scans[scan.uuid] = scan

    plot_received = []
    EventBroker().register(PlotFocusEvent, plot_received.append)
    EventBroker().publish(FocusEvent(ids=[scan.uuid]))

    assert len(plot_received) == 0


# ---------------------------------------------------------------------------
# _handle_focus_event — plot routing
# ---------------------------------------------------------------------------


def test_handle_focus_event_plot_publishes_plot_focus_event(model):
    plot = make_plot()
    scan = make_raw_scan()
    model.tavi_data.raw_scans[scan.uuid] = scan
    model.tavi_data.plots[plot.uuid] = plot

    received = []
    EventBroker().register(PlotFocusEvent, received.append)
    EventBroker().publish(FocusEvent(ids=[plot.uuid]))

    assert len(received) == 1
    assert received[0].plots[0].uuid == plot.uuid


def test_handle_focus_event_unmatched_uuid_is_noop_not_a_crash(model):
    """An unsaved preview plot's uuid belongs to the plot model, never to TaviData — must not raise."""
    received_raw = []
    received_plot = []
    EventBroker().register(RawScanFocusEvent, received_raw.append)
    EventBroker().register(PlotFocusEvent, received_plot.append)

    EventBroker().publish(FocusEvent(ids=[UUID(value="not-persisted")]))

    assert received_raw == []
    assert received_plot == []


def test_handle_focus_event_plot_publishes_scans_referenced_by_its_series(model):
    plot = make_plot()
    scan = make_raw_scan()
    model.tavi_data.raw_scans[scan.uuid] = scan
    model.tavi_data.plots[plot.uuid] = plot

    received = []
    EventBroker().register(PlotFocusEvent, received.append)
    EventBroker().publish(FocusEvent(ids=[plot.uuid]))

    assert received[0].scans[plot.series[0].source_scan_uuid].uuid == scan.uuid


def test_handle_focus_event_plot_does_not_publish_raw_scan_event(model):
    plot = make_plot()
    scan = make_raw_scan()
    model.tavi_data.raw_scans[scan.uuid] = scan
    model.tavi_data.plots[plot.uuid] = plot

    raw_received = []
    EventBroker().register(RawScanFocusEvent, raw_received.append)
    EventBroker().publish(FocusEvent(ids=[plot.uuid]))

    assert len(raw_received) == 0


# ---------------------------------------------------------------------------
# _handle_focus_event — mixed routing
# ---------------------------------------------------------------------------


def test_handle_focus_event_mixed_uuids_publishes_both_events(model):
    scan = make_raw_scan()
    plot = make_plot()
    model.tavi_data.raw_scans[scan.uuid] = scan
    model.tavi_data.plots[plot.uuid] = plot

    raw_received = []
    plot_received = []
    EventBroker().register(RawScanFocusEvent, raw_received.append)
    EventBroker().register(PlotFocusEvent, plot_received.append)

    EventBroker().publish(FocusEvent(ids=[scan.uuid, plot.uuid]))

    assert len(raw_received) == 1
    assert len(plot_received) == 1


def test_handle_focus_event_empty_ids_publishes_nothing(model):
    raw_received = []
    plot_received = []
    EventBroker().register(RawScanFocusEvent, raw_received.append)
    EventBroker().register(PlotFocusEvent, plot_received.append)

    EventBroker().publish(FocusEvent(ids=[]))

    assert len(raw_received) == 0
    assert len(plot_received) == 0


# ---------------------------------------------------------------------------
# _handle_active_plot_focus_event
# ---------------------------------------------------------------------------


def test_init_registers_active_plot_focus_event_handler(model):
    broker = EventBroker()
    assert model._handle_active_plot_focus_event in broker.registry[FocusActivePlotEvent]


def test_handle_active_plot_focus_event_announces_matching_saved_plot(model):
    plot = make_plot()
    scan = make_raw_scan()
    model.tavi_data.raw_scans[scan.uuid] = scan
    model.tavi_data.plots[plot.uuid] = plot

    received = []
    EventBroker().register(ActivePlotChangedEvent, received.append)
    EventBroker().publish(FocusActivePlotEvent(uuid=plot.uuid))

    assert len(received) == 1
    assert received[0].plot.uuid == plot.uuid


def test_handle_active_plot_focus_event_does_not_republish_plot_focus_event(model):
    """Switching the active plot must not re-resolve or re-render the whole focused batch."""
    plot = make_plot()
    scan = make_raw_scan()
    model.tavi_data.raw_scans[scan.uuid] = scan
    model.tavi_data.plots[plot.uuid] = plot

    received = []
    EventBroker().register(PlotFocusEvent, received.append)
    EventBroker().publish(FocusActivePlotEvent(uuid=plot.uuid))

    assert received == []


def test_handle_active_plot_focus_event_ignores_unmatched_uuid(model):
    """A uuid this model doesn't own (e.g. an unsaved preview plot's) must be a no-op here."""
    received = []
    EventBroker().register(ActivePlotChangedEvent, received.append)

    EventBroker().publish(FocusActivePlotEvent(uuid=UUID(value="not-persisted")))

    assert received == []


def test_handle_active_plot_focus_event_ignores_raw_scan_uuid(model):
    """A raw-scan uuid resolves to a RawScan, not a Plot — must not be announced as an active plot."""
    scan = make_raw_scan()
    model.tavi_data.raw_scans[scan.uuid] = scan

    received = []
    EventBroker().register(ActivePlotChangedEvent, received.append)

    EventBroker().publish(FocusActivePlotEvent(uuid=scan.uuid))

    assert received == []


# ---------------------------------------------------------------------------
# _handle_save_plot_event
# ---------------------------------------------------------------------------


def test_init_registers_save_plot_event_handler(model):
    broker = EventBroker()
    assert model._handle_save_plot_event in broker.registry[SavePlotEvent]


def test_handle_save_plot_event_stores_plot(model):
    plot = make_plot()

    EventBroker().publish(SavePlotEvent(plot=plot))

    assert model.tavi_data.plots[plot.uuid].uuid == plot.uuid


def test_handle_save_plot_event_publishes_plot_append_event(model):
    plot = make_plot()

    received = []
    EventBroker().register(PlotAppendEvent, received.append)
    EventBroker().publish(SavePlotEvent(plot=plot))

    assert len(received) == 1
    assert received[0].uuid == plot.uuid
    assert received[0].friendly_path == ""


def test_handle_save_plot_event_friendly_name_is_run_name_plus_plot_suffix(model):
    plot = make_plot(series=[make_series("run1")])

    received = []
    EventBroker().register(PlotAppendEvent, received.append)
    EventBroker().publish(SavePlotEvent(plot=plot))

    assert received[0].friendly_name == "run1_Plot"


def test_handle_save_plot_event_friendly_name_concatenates_multiple_run_names(model):
    plot = make_plot(series=[make_series("run1"), make_series("run2")])

    received = []
    EventBroker().register(PlotAppendEvent, received.append)
    EventBroker().publish(SavePlotEvent(plot=plot))

    assert received[0].friendly_name == "run1_run2_Plot"
