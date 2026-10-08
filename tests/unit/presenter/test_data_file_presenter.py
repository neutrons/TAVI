"""Tests for DataFilePresenter."""

from unittest.mock import MagicMock

import pytest

from tavi.frontend.presenter.data_file_presenter import DataFilePresenter
from tavi.frontend.view.data_file_view import DataFileView
from tavi.library.data.plot import PlotSeries
from tavi.library.data.scan import UUID, Provenance, RawScan, ScanData, ScanMetadata, TaviMetadata
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.presenter_event import ClearFocusEvent, SyncStageEvent


def make_scan(uuid_val="scan-001", data=None, metadata=None) -> RawScan:
    if data is None:
        data = {"qh": [1.0, 2.0], "en": [3.0, 4.0]}
    if metadata is None:
        metadata = ScanMetadata(
            data={"scan": "1"},
            categories={"ORNL Metadata": ["scan"]},
        )
    return RawScan(
        uuid=UUID(value=uuid_val),
        data=ScanData(data=data),
        metadata=metadata,
        tavimeta=TaviMetadata(default_axis=("qh", "en"), friendly_name="test_scan", friendly_path="/exp1"),
        prov=Provenance(raw_file="scan.dat", contributing_scans={UUID(value=uuid_val): 1}),
    )


def make_series(uuid_val="scan-001") -> PlotSeries:
    return PlotSeries(
        source_scan_uuid=UUID(value=uuid_val),
        scan_name="test_scan",
        normalized_by=None,
        x_name="qh",
        y_name="en",
        error_name="err",
    )


def make_stage_event(*scans: RawScan) -> SyncStageEvent:
    """A ``SyncStageEvent`` staging one series per scan, in order - the first leads."""
    return SyncStageEvent(
        series=[make_series(scan.uuid.value) for scan in scans], scans={scan.uuid: scan for scan in scans}
    )


@pytest.fixture
def presenter(qtbot):
    p = DataFilePresenter()
    qtbot.addWidget(p._view)
    return p


def test_init_view_is_data_file_view(presenter):
    assert isinstance(presenter._view, DataFileView)


def test_init_registers_clear_focus_and_sync_stage(presenter):
    broker = EventBroker()
    assert presenter.handle_clear_focus in broker.registry[ClearFocusEvent]
    assert presenter.handle_sync_stage in broker.registry[SyncStageEvent]


# ---------------------------------------------------------------------------
# handle_sync_stage
# ---------------------------------------------------------------------------


def test_sync_stage_populates_columns_from_the_lead_scan(presenter):
    scan = make_scan(data={"qh": [1.0, 2.0], "en": [3.0, 4.0]})
    presenter._view.populate_columns = MagicMock()

    presenter.handle_sync_stage(make_stage_event(scan))

    presenter._view.populate_columns.assert_called_once_with(scan.data.data)


def test_sync_stage_populates_variables_with_column_names(presenter):
    scan = make_scan(data={"qh": [1.0], "en": [2.0]})
    presenter._view.populate_variables = MagicMock()

    presenter.handle_sync_stage(make_stage_event(scan))

    args = presenter._view.populate_variables.call_args.args
    assert set(args[0]) == {"qh", "en"}


def test_sync_stage_populates_metadata_by_category(presenter):
    metadata = ScanMetadata(
        data={"scan": "1", "proposal": "9865"},
        categories={"ORNL Metadata": ["scan", "proposal"]},
    )
    scan = make_scan(metadata=metadata)
    presenter._view.populate_metadata = MagicMock()

    presenter.handle_sync_stage(make_stage_event(scan))

    presenter._view.populate_metadata.assert_called_once_with(metadata.by_category())


def test_sync_stage_shows_only_the_lead_scan(presenter):
    scan1 = make_scan(uuid_val="scan-001", data={"qh": [1.0]})
    scan2 = make_scan(uuid_val="scan-002", data={"qh": [2.0]})
    presenter._view.populate_columns = MagicMock()

    presenter.handle_sync_stage(make_stage_event(scan1, scan2))

    presenter._view.populate_columns.assert_called_once_with(scan1.data.data)


def test_sync_stage_sets_title_from_friendly_name(presenter):
    received = []
    presenter._view.title_changed.connect(received.append)

    presenter.handle_sync_stage(make_stage_event(make_scan()))

    assert received == ["Data File (test_scan)"]


def test_sync_stage_with_nothing_staged_clears_view(presenter):
    presenter._view.clear_data = MagicMock()

    presenter.handle_sync_stage(SyncStageEvent(series=[]))

    presenter._view.clear_data.assert_called_once()


def test_sync_stage_via_event_broker(presenter):
    presenter._view.populate_columns = MagicMock()

    EventBroker().publish(make_stage_event(make_scan()))

    presenter._view.populate_columns.assert_called_once()


# ---------------------------------------------------------------------------
# handle_clear_focus
# ---------------------------------------------------------------------------


def test_clear_focus_clears_view_and_resets_title(presenter):
    presenter._view.clear_data = MagicMock()
    received = []
    presenter._view.title_changed.connect(received.append)

    EventBroker().publish(ClearFocusEvent())

    presenter._view.clear_data.assert_called_once()
    assert received == ["Data File"]
