"""Tests for CombinePresenter."""

from unittest.mock import MagicMock

import pytest

from tavi.frontend.presenter.combine_presenter import CombinePresenter
from tavi.frontend.view.combine_view import CombineView
from tavi.library.data.enum.rebin_mode import RebinMode
from tavi.library.data.scan import UUID, Provenance, RawScan, ScanData, ScanMetadata, TaviMetadata
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.model_event import RemoveRawScanEvent
from tavi.meta.event.type.presenter_event import ClearFocusEvent, FocusRawScanEvent


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


def focus(*scans: RawScan) -> None:
    EventBroker().publish(FocusRawScanEvent(scans=list(scans)))


@pytest.fixture
def model():
    return MagicMock()


@pytest.fixture
def presenter(qtbot, model):
    p = CombinePresenter(model)
    qtbot.addWidget(p._view)
    return p


def last_request(model):
    """The CombineRequest the model was last asked to run."""
    return model.combine.call_args.args[0]


# ---------------------------------------------------------------------------
# construction
# ---------------------------------------------------------------------------


def test_init_view_is_combine_view(presenter):
    assert isinstance(presenter._view, CombineView)


def test_nothing_is_staged_initially(presenter):
    assert presenter._view.scan_names() == []


# ---------------------------------------------------------------------------
# staging focused scans
# ---------------------------------------------------------------------------


def test_add_scans_lists_the_focused_scans(presenter):
    focus(make_scan("a", "scan0001"), make_scan("b", "scan0002"))

    presenter._view.add_scans_btn.click()

    assert presenter._view.scan_names() == ["scan0001", "scan0002"]


def test_add_scans_with_nothing_focused_lists_nothing(presenter):
    presenter._view.add_scans_btn.click()

    assert presenter._view.scan_names() == []


def test_add_scans_keeps_what_an_earlier_selection_staged(presenter):
    focus(make_scan("a", "scan0001"))
    presenter._view.add_scans_btn.click()

    EventBroker().publish(ClearFocusEvent())
    focus(make_scan("b", "scan0002"))
    presenter._view.add_scans_btn.click()

    assert presenter._view.scan_names() == ["scan0001", "scan0002"]


def test_clear_focus_leaves_the_staged_list_alone(presenter):
    focus(make_scan("a", "scan0001"))
    presenter._view.add_scans_btn.click()

    EventBroker().publish(ClearFocusEvent())

    assert presenter._view.scan_names() == ["scan0001"]


def test_add_scans_after_clear_focus_stages_nothing_again(presenter):
    focus(make_scan("a", "scan0001"))
    EventBroker().publish(ClearFocusEvent())

    presenter._view.add_scans_btn.click()

    assert presenter._view.scan_names() == []


def test_focus_offers_the_scan_columns_in_the_column_fields(presenter):
    focus(make_scan("a", "scan0001"))

    items = [presenter._view.x_axis_combo.itemText(i) for i in range(presenter._view.x_axis_combo.count())]
    assert items == ["qh", "detector", "mcu"]


def test_focus_seeds_the_axis_fields_from_the_scan_default_axis(presenter):
    """Otherwise the panel keeps naming columns the newly focused scans may not carry."""
    focus(make_scan("a", "scan0001"))

    assert (presenter._view.x_axis(), presenter._view.y_axis()) == ("qh", "detector")


def test_focus_leaves_the_normalization_at_mcu(presenter):
    focus(make_scan("a", "scan0001"))

    assert presenter._view.normalization() == "mcu"


def test_a_later_focus_reseeds_the_axis_fields(presenter):
    focus(make_scan("a", "scan0001"))
    other = make_scan("b", "scan0002", columns={"en": [1.0], "det_2": [2.0]})
    other.tavimeta.default_axis = ("en", "det_2")

    focus(other)

    assert (presenter._view.x_axis(), presenter._view.y_axis()) == ("en", "det_2")


# ---------------------------------------------------------------------------
# removing
# ---------------------------------------------------------------------------


def test_removing_a_scan_drops_it_from_the_next_request(presenter, model):
    focus(make_scan("a", "scan0001"), make_scan("b", "scan0002"))
    presenter._view.add_scans_btn.click()
    presenter._view.scan_list.item(0).setSelected(True)
    presenter._view.remove_scans_btn.click()

    presenter._view.append_btn.click()

    assert last_request(model).source_scan_uuids == [UUID(value="b")]


def test_a_removed_raw_scan_leaves_the_list(presenter):
    focus(make_scan("a", "scan0001"), make_scan("b", "scan0002"))
    presenter._view.add_scans_btn.click()

    EventBroker().publish(RemoveRawScanEvent(uuid=UUID(value="a")))

    assert presenter._view.scan_names() == ["scan0002"]


def test_a_removed_raw_scan_drops_out_of_the_next_request(presenter, model):
    focus(make_scan("a", "scan0001"), make_scan("b", "scan0002"))
    presenter._view.add_scans_btn.click()
    EventBroker().publish(RemoveRawScanEvent(uuid=UUID(value="a")))

    presenter._view.append_btn.click()

    assert last_request(model).source_scan_uuids == [UUID(value="b")]


def test_a_removed_raw_scan_cannot_be_staged_again(presenter):
    focus(make_scan("a", "scan0001"))
    EventBroker().publish(RemoveRawScanEvent(uuid=UUID(value="a")))

    presenter._view.add_scans_btn.click()

    assert presenter._view.scan_names() == []


# ---------------------------------------------------------------------------
# Append
# ---------------------------------------------------------------------------


def test_append_asks_the_model_to_combine(presenter, model):
    focus(make_scan("a", "scan0001"))
    presenter._view.add_scans_btn.click()

    presenter._view.append_btn.click()

    model.combine.assert_called_once()


def test_append_sends_no_rebin(presenter, model):
    focus(make_scan("a", "scan0001"))
    presenter._view.add_scans_btn.click()

    presenter._view.append_btn.click()

    assert last_request(model).rebin_mode is RebinMode.NONE


def test_append_sends_no_normalization_or_scale(presenter, model):
    """An append applies no weighting, and CombineOps refuses one it cannot honour."""
    focus(make_scan("a", "scan0001"))
    presenter._view.add_scans_btn.click()
    presenter._view.normalization_combo.setEditText("mcu")
    presenter._view.multiplication_edit.setText("10")

    presenter._view.append_btn.click()

    request = last_request(model)
    assert request.normalization == ""
    assert request.multiplication == "1"


def test_append_sends_the_axis_fields(presenter, model):
    focus(make_scan("a", "scan0001"))
    presenter._view.add_scans_btn.click()
    presenter._view.x_axis_combo.setEditText("qh")
    presenter._view.y_axis_combo.setEditText("detector")

    presenter._view.append_btn.click()

    request = last_request(model)
    assert (request.x_axis, request.y_axis) == ("qh", "detector")


# ---------------------------------------------------------------------------
# Combine
# ---------------------------------------------------------------------------


def test_combine_sends_the_equal_step_tab(presenter, model):
    focus(make_scan("a", "scan0001"))
    presenter._view.add_scans_btn.click()
    presenter._view.rebin_tabs.setCurrentIndex(0)
    presenter._view.rebin_start_edit.setText("-1")
    presenter._view.rebin_stop_edit.setText("3")
    presenter._view.rebin_step_edit.setText("0.1")

    presenter._view.combine_btn.click()

    request = last_request(model)
    assert request.rebin_mode is RebinMode.EQUAL_STEP
    assert (request.start, request.finish, request.step) == ("-1", "3", "0.1")


def test_combine_sends_the_tolerance_tab(presenter, model):
    focus(make_scan("a", "scan0001"))
    presenter._view.add_scans_btn.click()
    presenter._view.rebin_tabs.setCurrentIndex(1)
    presenter._view.rebin_tolerance_edit.setText("0.005")

    presenter._view.combine_btn.click()

    request = last_request(model)
    assert request.rebin_mode is RebinMode.TOLERANCE
    assert request.tolerance == "0.005"


def test_combine_sends_the_normalization_and_scale(presenter, model):
    focus(make_scan("a", "scan0001"))
    presenter._view.add_scans_btn.click()
    presenter._view.multiplication_edit.setText("10000")

    presenter._view.combine_btn.click()

    request = last_request(model)
    assert request.normalization == "mcu"
    assert request.multiplication == "10000"


def test_combine_sends_the_staged_scans_in_order(presenter, model):
    focus(make_scan("a", "scan0001"), make_scan("b", "scan0002"))
    presenter._view.add_scans_btn.click()

    presenter._view.combine_btn.click()

    assert last_request(model).source_scan_uuids == [UUID(value="a"), UUID(value="b")]
