import numpy as np
import pytest
import unittest
from unittest import mock

from tavi.backend.model.plot_model import PlotModel
from tavi.library.data.plot import Plot
from tavi.library.data.scan import UUID, RawScan, ScanData, ScanMetadata, TaviMetadata, Provenance
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.presenter_event import PlotFocusEvent, RawScanFocusEvent


def make_raw_scan(x_col="qh", x_vals=None, y_col="en", y_vals=None, norm=("monitor", 1.0)) -> RawScan:
    if x_vals is None:
        x_vals = [1.0, 2.0, 3.0]
    if y_vals is None:
        y_vals = [1.0, 2.0, 3.0]
    return RawScan(
        uuid=UUID(value="scan-001"),
        data=ScanData(data={x_col: x_vals, y_col: y_vals}),
        metadata=ScanMetadata(),
        tavimeta=TaviMetadata(
            default_axis=(x_col, y_col),
            normalization=norm,
            friendly_name="test_scan",
            friendly_path="/test_path",
        ),
        prov=Provenance(raw_file="scan0001.dat", contributing_scans={UUID(value="scan-001"): 1}),
    )


class TestPlotModel(unittest.TestCase):
    def setUp(self):
        self.broker = EventBroker()
        self.model = PlotModel(plots=[])

    def tearDown(self):
        pass

    def test_registers_raw_scan_focus_event_handler(self):
        assert RawScanFocusEvent in self.broker.registry
        assert len(self.broker.registry[RawScanFocusEvent]) == 1

    def test_raw_scan_focus_event_publishes_plot_focus_event(self):
        received_events: list[PlotFocusEvent] = []

        self.broker.register(PlotFocusEvent, received_events.append)
        scan = make_raw_scan()
        self.broker.publish(RawScanFocusEvent(scans=[scan]))

        assert len(received_events) == 1
        assert isinstance(received_events[0], PlotFocusEvent)

    def test_plot_contains_correct_scan_name(self):
        received: list[PlotFocusEvent] = []
        self.broker.register(PlotFocusEvent, received.append)

        scan = make_raw_scan()
        self.broker.publish(RawScanFocusEvent(scans=[scan]))

        assert received[0].plots[0].scan_name == "test_scan"

    def test_plot_contains_correct_x_data(self):
        received: list[PlotFocusEvent] = []
        self.broker.register(PlotFocusEvent, received.append)

        scan = make_raw_scan(x_col="qh", x_vals=[10.0, 20.0, 30.0])
        self.broker.publish(RawScanFocusEvent(scans=[scan]))

        np.testing.assert_array_equal(received[0].plots[0].x, [10.0, 20.0, 30.0])
        assert received[0].plots[0].x_name == "qh"

    def test_plot_normalized_by_set_from_tavimeta(self):
        received: list[PlotFocusEvent] = []
        self.broker.register(PlotFocusEvent, received.append)

        scan = make_raw_scan(norm=("detector", 1.0))
        self.broker.publish(RawScanFocusEvent(scans=[scan]))

        assert received[0].plots[0].normalized_by == "detector"

    def test_plot_normalized_by_none_when_no_normalization(self):
        received: list[PlotFocusEvent] = []
        self.broker.register(PlotFocusEvent, received.append)

        scan = make_raw_scan(norm=None)
        self.broker.publish(RawScanFocusEvent(scans=[scan]))

        assert received[0].plots[0].normalized_by is None

    def test_plot_err_is_zeros(self):
        received: list[PlotFocusEvent] = []
        self.broker.register(PlotFocusEvent, received.append)

        scan = make_raw_scan(x_vals=[1.0, 2.0, 3.0])
        self.broker.publish(RawScanFocusEvent(scans=[scan]))

        np.testing.assert_array_equal(received[0].plots[0].err, [0.0, 0.0, 0.0])

    def test_only_first_scan_processed(self):
        """PlotModel only handles scans[0] from the event."""
        received: list[PlotFocusEvent] = []
        self.broker.register(PlotFocusEvent, received.append)

        scan1 = make_raw_scan(x_col="qh", x_vals=[1.0])
        scan2 = make_raw_scan(x_col="qh", x_vals=[99.0])
        # scan2 has different friendly_name setup — both have "test_scan"
        self.broker.publish(RawScanFocusEvent(scans=[scan1, scan2]))

        assert len(received[0].plots) == 1
