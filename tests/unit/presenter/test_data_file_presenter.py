"""Tests for DataFilePresenter."""

from unittest.mock import MagicMock

import pytest

from tavi.frontend.presenter.data_file_presenter import DataFilePresenter
from tavi.frontend.view.data_file_view import DataFileView
from tavi.library.data.scan import UUID, Provenance, RawScan, ScanData, ScanMetadata, TaviMetadata
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.presenter_event import RawScanFocusEvent


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


@pytest.fixture
def presenter(qtbot):
    p = DataFilePresenter()
    qtbot.addWidget(p._view)
    return p


def test_init_view_is_data_file_view(presenter):
    assert isinstance(presenter._view, DataFileView)


def test_init_registers_raw_scan_focus_event(presenter):
    broker = EventBroker()
    assert presenter.handle_raw_scan_focus in broker.registry[RawScanFocusEvent]


def test_handle_raw_scan_focus_no_scans_clears_view(presenter):
    presenter._view.clear_data = MagicMock()
    presenter._view.populate_columns = MagicMock()

    presenter.handle_raw_scan_focus(RawScanFocusEvent(scans=[]))

    presenter._view.clear_data.assert_called_once()
    presenter._view.populate_columns.assert_not_called()


def test_handle_raw_scan_focus_populates_columns_from_first_scan(presenter):
    scan = make_scan(data={"qh": [1.0, 2.0], "en": [3.0, 4.0]})
    presenter._view.populate_columns = MagicMock()

    presenter.handle_raw_scan_focus(RawScanFocusEvent(scans=[scan]))

    presenter._view.populate_columns.assert_called_once_with(scan.data.data)


def test_handle_raw_scan_focus_populates_variables_with_column_names(presenter):
    scan = make_scan(data={"qh": [1.0], "en": [2.0]})
    presenter._view.populate_variables = MagicMock()

    presenter.handle_raw_scan_focus(RawScanFocusEvent(scans=[scan]))

    args = presenter._view.populate_variables.call_args.args
    assert set(args[0]) == {"qh", "en"}


def test_handle_raw_scan_focus_populates_metadata_by_category(presenter):
    metadata = ScanMetadata(
        data={"scan": "1", "proposal": "9865"},
        categories={"ORNL Metadata": ["scan", "proposal"]},
    )
    scan = make_scan(metadata=metadata)
    presenter._view.populate_metadata = MagicMock()

    presenter.handle_raw_scan_focus(RawScanFocusEvent(scans=[scan]))

    presenter._view.populate_metadata.assert_called_once_with(metadata.by_category())


def test_handle_raw_scan_focus_uses_only_first_scan(presenter):
    scan1 = make_scan(uuid_val="scan-001", data={"qh": [1.0]})
    scan2 = make_scan(uuid_val="scan-002", data={"qh": [2.0]})
    presenter._view.populate_columns = MagicMock()

    presenter.handle_raw_scan_focus(RawScanFocusEvent(scans=[scan1, scan2]))

    presenter._view.populate_columns.assert_called_once_with(scan1.data.data)


def test_handle_raw_scan_focus_via_event_broker(presenter):
    scan = make_scan()
    presenter._view.populate_columns = MagicMock()

    EventBroker().publish(RawScanFocusEvent(scans=[scan]))

    presenter._view.populate_columns.assert_called_once_with(scan.data.data)
