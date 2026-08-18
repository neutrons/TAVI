import pytest
import unittest
from unittest import mock

from tavi.backend.model.plot_model import PlotModel
from tavi.library.data.plot import Plot, PlotFields, PlotSeries
from tavi.library.data.scan import UUID, RawScan, ScanData, ScanMetadata, TaviMetadata, Provenance
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.model_event import PlotAppendEvent
from tavi.meta.event.type.presenter_event import PlotFocusEvent, RawScanFocusEvent, SavePlotEvent


def make_plot_fields(**overrides) -> PlotFields:
    defaults = dict(
        y_axis="",
        x_axis="",
        rebin_mode="none",
        rebin_start="0",
        rebin_stop="2",
        rebin_step="0.02",
        preset_type="none",
        preset_channel="",
        preset_value="1",
    )
    defaults.update(overrides)
    return PlotFields(**defaults)


def make_raw_scan(x_col="qh", x_vals=None, y_col="en", y_vals=None, norm=("monitor", 1.0), uuid_val="scan-001") -> RawScan:
    if x_vals is None:
        x_vals = [1.0, 2.0, 3.0]
    if y_vals is None:
        y_vals = [1.0, 2.0, 3.0]
    return RawScan(
        uuid=UUID(value=uuid_val),
        data=ScanData(data={x_col: x_vals, y_col: y_vals}),
        metadata=ScanMetadata(),
        tavimeta=TaviMetadata(
            default_axis=(x_col, y_col),
            normalization=norm,
            friendly_name="test_scan",
            friendly_path="/test_path",
        ),
        prov=Provenance(raw_file="scan0001.dat", contributing_scans={UUID(value=uuid_val): 1}),
    )


def make_series(scan_name, uuid_val="scan-001") -> PlotSeries:
    return PlotSeries(
        source_scan_uuid=UUID(value=uuid_val),
        scan_name=scan_name,
        normalized_by=None,
        x_name="qh",
        y_name="en",
        error_name="err",
    )


def make_plot(uuid_val="plot-001", series=None) -> Plot:
    if series is None:
        series = [make_series("test_plot")]
    return Plot(uuid=UUID(value=uuid_val), series=series)


class TestPlotModel(unittest.TestCase):
    def setUp(self):
        self.broker = EventBroker()
        self.raw_scans = {}
        self.plots = {}
        self.model = PlotModel(plots=self.plots, raw_scans=self.raw_scans)

    def tearDown(self):
        pass

    def test_registers_raw_scan_focus_event_handler(self):
        assert RawScanFocusEvent in self.broker.registry
        assert len(self.broker.registry[RawScanFocusEvent]) == 1

    def test_raw_scan_focus_event_empty_scans_is_noop(self):
        """No scan selected: must not IndexError, and must not publish a PlotFocusEvent."""
        received_events: list[PlotFocusEvent] = []
        self.broker.register(PlotFocusEvent, received_events.append)

        self.broker.publish(RawScanFocusEvent(scans=[]))

        assert len(received_events) == 0
        assert self.model._last_plot is None

    def test_raw_scan_focus_event_publishes_plot_focus_event(self):
        received_events: list[PlotFocusEvent] = []

        self.broker.register(PlotFocusEvent, received_events.append)
        scan = make_raw_scan()
        self.raw_scans[scan.uuid] = scan
        self.broker.publish(RawScanFocusEvent(scans=[scan]))

        assert len(received_events) == 1
        assert isinstance(received_events[0], PlotFocusEvent)

    def test_plot_has_single_series_with_correct_scan_name(self):
        received: list[PlotFocusEvent] = []
        self.broker.register(PlotFocusEvent, received.append)

        scan = make_raw_scan()
        self.raw_scans[scan.uuid] = scan
        self.broker.publish(RawScanFocusEvent(scans=[scan]))

        assert len(received[0].plots[0].series) == 1
        assert received[0].plots[0].series[0].scan_name == "test_scan"

    def test_plot_series_points_at_source_scan(self):
        received: list[PlotFocusEvent] = []
        self.broker.register(PlotFocusEvent, received.append)

        scan = make_raw_scan(x_col="qh")
        self.raw_scans[scan.uuid] = scan
        self.broker.publish(RawScanFocusEvent(scans=[scan]))

        series = received[0].plots[0].series[0]
        assert series.source_scan_uuid == scan.uuid
        assert series.x_name == "qh"

    def test_plot_normalized_by_not_applied_by_default_even_when_tavimeta_has_normalization(self):
        received: list[PlotFocusEvent] = []
        self.broker.register(PlotFocusEvent, received.append)

        scan = make_raw_scan(norm=("detector", 1.0))
        self.raw_scans[scan.uuid] = scan
        self.broker.publish(RawScanFocusEvent(scans=[scan]))

        assert received[0].plots[0].series[0].normalized_by is None

    def test_plot_normalized_by_none_when_no_normalization(self):
        received: list[PlotFocusEvent] = []
        self.broker.register(PlotFocusEvent, received.append)

        scan = make_raw_scan(norm=None)
        self.raw_scans[scan.uuid] = scan
        self.broker.publish(RawScanFocusEvent(scans=[scan]))

        assert received[0].plots[0].series[0].normalized_by is None

    def test_only_first_scan_processed(self):
        """PlotModel only handles scans[0] from the event."""
        received: list[PlotFocusEvent] = []
        self.broker.register(PlotFocusEvent, received.append)

        scan1 = make_raw_scan(x_col="qh", x_vals=[1.0])
        scan2 = make_raw_scan(x_col="qh", x_vals=[99.0])
        # scan2 has different friendly_name setup — both have "test_scan"
        self.raw_scans[scan1.uuid] = scan1
        self.broker.publish(RawScanFocusEvent(scans=[scan1, scan2]))

        assert len(received[0].plots) == 1

    def test_plot_focus_event_carries_the_referenced_scan(self):
        received: list[PlotFocusEvent] = []
        self.broker.register(PlotFocusEvent, received.append)

        scan = make_raw_scan()
        self.raw_scans[scan.uuid] = scan
        self.broker.publish(RawScanFocusEvent(scans=[scan]))

        assert received[0].scans[scan.uuid].uuid == scan.uuid

    def test_update_fields_updates_axis_names_on_series(self):
        scan = make_raw_scan(x_col="qh", y_col="en")
        self.raw_scans[scan.uuid] = scan
        self.broker.register(RawScanFocusEvent, self.model._handle_raw_scan_focus_event)
        self.broker.publish(RawScanFocusEvent(scans=[scan]))

        received: list[PlotFocusEvent] = []
        self.broker.register(PlotFocusEvent, received.append)

        response = self.model.update_fields(make_plot_fields(x_axis="en", y_axis="qh"))

        assert response.code.name == "OK"
        series = received[0].plots[0].series[0]
        assert series.x_name == "en"
        assert series.y_name == "qh"

    def test_update_fields_stores_friendly_name_on_plot(self):
        scan = make_raw_scan(x_col="qh", y_col="en")
        self.raw_scans[scan.uuid] = scan
        self.model._handle_raw_scan_focus_event(RawScanFocusEvent(scans=[scan]))

        received: list[PlotFocusEvent] = []
        self.broker.register(PlotFocusEvent, received.append)

        self.model.update_fields(make_plot_fields(x_axis="en", y_axis="qh", friendly_name="my_custom_name"))

        assert received[0].plots[0].friendly_name == "my_custom_name"

    def test_update_fields_blank_friendly_name_leaves_it_unset(self):
        scan = make_raw_scan(x_col="qh", y_col="en")
        self.raw_scans[scan.uuid] = scan
        self.model._handle_raw_scan_focus_event(RawScanFocusEvent(scans=[scan]))

        received: list[PlotFocusEvent] = []
        self.broker.register(PlotFocusEvent, received.append)

        self.model.update_fields(make_plot_fields(x_axis="en", y_axis="qh", friendly_name="   "))

        assert received[0].plots[0].friendly_name is None

    def test_update_fields_no_focused_plot_is_noop(self):
        response = self.model.update_fields(make_plot_fields())

        assert response.code.name == "OK"

    def test_update_fields_unknown_column_is_noop(self):
        scan = make_raw_scan(x_col="qh", y_col="en")
        self.raw_scans[scan.uuid] = scan
        self.model._handle_raw_scan_focus_event(RawScanFocusEvent(scans=[scan]))

        received: list[PlotFocusEvent] = []
        self.broker.register(PlotFocusEvent, received.append)

        response = self.model.update_fields(make_plot_fields(x_axis="nonexistent", y_axis="en"))

        assert response.code.name == "OK"
        assert len(received) == 0

    def test_raw_scan_focus_normalized_by_value_not_applied_by_default(self):
        received: list[PlotFocusEvent] = []
        self.broker.register(PlotFocusEvent, received.append)

        scan = make_raw_scan(norm=("monitor", 2.5))
        self.raw_scans[scan.uuid] = scan
        self.broker.publish(RawScanFocusEvent(scans=[scan]))

        series = received[0].plots[0].series[0]
        assert series.normalized_by is None
        assert series.normalized_by_value is None

    def test_raw_scan_focus_normalized_by_value_none_when_no_normalization(self):
        received: list[PlotFocusEvent] = []
        self.broker.register(PlotFocusEvent, received.append)

        scan = make_raw_scan(norm=None)
        self.raw_scans[scan.uuid] = scan
        self.broker.publish(RawScanFocusEvent(scans=[scan]))

        assert received[0].plots[0].series[0].normalized_by_value is None

    def test_update_fields_normalize_sets_channel_and_value(self):
        scan = make_raw_scan(x_col="qh", y_col="en", norm=None)
        self.raw_scans[scan.uuid] = scan
        self.model._handle_raw_scan_focus_event(RawScanFocusEvent(scans=[scan]))

        received: list[PlotFocusEvent] = []
        self.broker.register(PlotFocusEvent, received.append)

        response = self.model.update_fields(
            make_plot_fields(x_axis="qh", y_axis="en", preset_type="normalize", preset_channel="qh", preset_value="3.5")
        )

        assert response.code.name == "OK"
        series = received[0].plots[0].series[0]
        assert series.normalized_by == "qh"
        assert series.normalized_by_value == 3.5

    def test_update_fields_normalize_unknown_channel_is_noop(self):
        scan = make_raw_scan(x_col="qh", y_col="en", norm=None)
        self.raw_scans[scan.uuid] = scan
        self.model._handle_raw_scan_focus_event(RawScanFocusEvent(scans=[scan]))

        received: list[PlotFocusEvent] = []
        self.broker.register(PlotFocusEvent, received.append)

        response = self.model.update_fields(
            make_plot_fields(x_axis="qh", y_axis="en", preset_type="normalize", preset_channel="nonexistent", preset_value="1")
        )

        assert response.code.name == "OK"
        assert len(received) == 0

    def test_update_fields_normalize_invalid_value_is_noop(self):
        scan = make_raw_scan(x_col="qh", y_col="en", norm=None)
        self.raw_scans[scan.uuid] = scan
        self.model._handle_raw_scan_focus_event(RawScanFocusEvent(scans=[scan]))

        received: list[PlotFocusEvent] = []
        self.broker.register(PlotFocusEvent, received.append)

        response = self.model.update_fields(
            make_plot_fields(x_axis="qh", y_axis="en", preset_type="normalize", preset_channel="qh", preset_value="not_a_number")
        )

        assert response.code.name == "OK"
        assert len(received) == 0

    def test_update_fields_none_preset_type_clears_normalization(self):
        scan = make_raw_scan(x_col="qh", y_col="en", norm=("monitor", 1.0))
        self.raw_scans[scan.uuid] = scan
        self.model._handle_raw_scan_focus_event(RawScanFocusEvent(scans=[scan]))

        received: list[PlotFocusEvent] = []
        self.broker.register(PlotFocusEvent, received.append)

        response = self.model.update_fields(make_plot_fields(x_axis="qh", y_axis="en", preset_type="none"))

        assert response.code.name == "OK"
        series = received[0].plots[0].series[0]
        assert series.normalized_by is None
        assert series.normalized_by_value is None

    def test_init_registers_save_plot_event_handler(self):
        assert self.model._handle_save_plot_event in self.broker.registry[SavePlotEvent]

    def test_handle_save_plot_event_stores_plot(self):
        plot = make_plot()

        self.broker.publish(SavePlotEvent(plot=plot))

        assert self.plots[plot.uuid].uuid == plot.uuid

    def test_handle_save_plot_event_publishes_plot_append_event(self):
        plot = make_plot()

        received = []
        self.broker.register(PlotAppendEvent, received.append)
        self.broker.publish(SavePlotEvent(plot=plot))

        assert len(received) == 1
        assert received[0].uuid == plot.uuid
        assert received[0].friendly_path == ""

    def test_handle_save_plot_event_friendly_name_is_run_name_plus_plot_suffix(self):
        plot = make_plot(series=[make_series("run1")])

        received = []
        self.broker.register(PlotAppendEvent, received.append)
        self.broker.publish(SavePlotEvent(plot=plot))

        assert received[0].friendly_name == "run1_Plot"

    def test_handle_save_plot_event_friendly_name_concatenates_multiple_run_names(self):
        plot = make_plot(series=[make_series("run1"), make_series("run2")])

        received = []
        self.broker.register(PlotAppendEvent, received.append)
        self.broker.publish(SavePlotEvent(plot=plot))

        assert received[0].friendly_name == "run1_run2_Plot"

    def test_handle_save_plot_event_second_save_of_same_run_name_increments_suffix(self):
        first = make_plot(uuid_val="plot-001", series=[make_series("run1")])
        second = make_plot(uuid_val="plot-002", series=[make_series("run1")])

        self.broker.publish(SavePlotEvent(plot=first))

        received = []
        self.broker.register(PlotAppendEvent, received.append)
        self.broker.publish(SavePlotEvent(plot=second))

        assert received[0].friendly_name == "run1_Plot(1)"

    def test_handle_save_plot_event_stamps_final_name_onto_stored_plot(self):
        first = make_plot(uuid_val="plot-001", series=[make_series("run1")])
        second = make_plot(uuid_val="plot-002", series=[make_series("run1")])

        self.broker.publish(SavePlotEvent(plot=first))
        self.broker.publish(SavePlotEvent(plot=second))

        assert self.plots[first.uuid].friendly_name == "run1_Plot"
        assert self.plots[second.uuid].friendly_name == "run1_Plot(1)"

    def test_handle_save_plot_event_user_provided_name_overrides_run_name_default(self):
        plot = make_plot(series=[make_series("run1")])
        plot = plot.model_copy(update={"friendly_name": "my_custom_name"})

        received = []
        self.broker.register(PlotAppendEvent, received.append)
        self.broker.publish(SavePlotEvent(plot=plot))

        assert received[0].friendly_name == "my_custom_name"

    def test_handle_save_plot_event_user_provided_name_still_deduped_against_existing(self):
        first = make_plot(uuid_val="plot-001", series=[make_series("run1")])
        first = first.model_copy(update={"friendly_name": "my_custom_name"})
        second = make_plot(uuid_val="plot-002", series=[make_series("run2")])
        second = second.model_copy(update={"friendly_name": "my_custom_name"})

        self.broker.publish(SavePlotEvent(plot=first))

        received = []
        self.broker.register(PlotAppendEvent, received.append)
        self.broker.publish(SavePlotEvent(plot=second))

        assert received[0].friendly_name == "my_custom_name(1)"
