"""Tests for tavi.frontend.view.combine_view."""

import pytest
from qtpy.QtWidgets import QGroupBox, QHBoxLayout, QLabel, QListWidget, QVBoxLayout

from tavi.frontend.view.combine_view import CombineView
from tavi.library.data.enum.rebin_mode import RebinMode


@pytest.fixture
def view(qtbot):
    w = CombineView()
    qtbot.addWidget(w)
    return w


def test_combine_view_starts_empty(view):
    assert view.scan_names() == []


def test_add_scans_lists_a_single_name(view):
    view.add_scans(["scan0001"])

    assert view.scan_names() == ["scan0001"]


def test_add_scans_lists_a_multiselection_in_order(view):
    view.add_scans(["scan0001", "scan0002", "scan0003"])

    assert view.scan_names() == ["scan0001", "scan0002", "scan0003"]


def test_add_scans_appends_to_what_is_already_listed(view):
    view.add_scans(["scan0001"])
    view.add_scans(["scan0002"])

    assert view.scan_names() == ["scan0001", "scan0002"]


def test_add_scans_skips_names_already_listed(view):
    view.add_scans(["scan0001", "scan0002"])
    view.add_scans(["scan0002", "scan0003"])

    assert view.scan_names() == ["scan0001", "scan0002", "scan0003"]


def test_add_scans_skips_duplicates_within_one_call(view):
    view.add_scans(["scan0001", "scan0001"])

    assert view.scan_names() == ["scan0001"]


def test_add_scans_with_nothing_selected_changes_nothing(view):
    view.add_scans(["scan0001"])
    view.add_scans([])

    assert view.scan_names() == ["scan0001"]


def test_clear_scans_empties_the_list(view):
    view.add_scans(["scan0001", "scan0002"])
    view.clear_scans()

    assert view.scan_names() == []


def test_add_scans_button_emits_add_scans_clicked(view, qtbot):
    with qtbot.waitSignal(view.add_scans_clicked, timeout=1000):
        view.add_scans_btn.click()


def test_add_scans_button_label(view):
    assert view.add_scans_btn.text() == "Add Scans"


# ---------------------------------------------------------------------------
# CombineView — removing scans
# ---------------------------------------------------------------------------


def _select(view, *names):
    for row in range(view.scan_list.count()):
        item = view.scan_list.item(row)
        item.setSelected(item.text() in names)


def test_remove_scans_button_label(view):
    assert view.remove_scans_btn.text() == "Remove Scans"


def test_remove_scans_button_disabled_while_nothing_is_selected(view):
    view.add_scans(["scan0001"])

    assert not view.remove_scans_btn.isEnabled()


def test_remove_scans_button_enabled_once_a_scan_is_selected(view):
    view.add_scans(["scan0001"])
    _select(view, "scan0001")

    assert view.remove_scans_btn.isEnabled()


def test_remove_scans_button_removes_a_single_selected_scan(view):
    view.add_scans(["scan0001", "scan0002"])
    _select(view, "scan0001")

    view.remove_scans_btn.click()

    assert view.scan_names() == ["scan0002"]


def test_remove_scans_button_removes_a_multiselection(view):
    view.add_scans(["scan0001", "scan0002", "scan0003"])
    _select(view, "scan0001", "scan0003")

    view.remove_scans_btn.click()

    assert view.scan_names() == ["scan0002"]


def test_remove_scans_button_emits_the_names_it_removed(view, qtbot):
    view.add_scans(["scan0001", "scan0002"])
    _select(view, "scan0002")

    with qtbot.waitSignal(view.remove_scans_clicked, timeout=1000) as blocker:
        view.remove_scans_btn.click()

    assert blocker.args == [["scan0002"]]


def test_remove_scans_button_disables_itself_again_after_removing(view):
    view.add_scans(["scan0001"])
    _select(view, "scan0001")

    view.remove_scans_btn.click()

    assert not view.remove_scans_btn.isEnabled()


def test_remove_scans_by_name(view):
    view.add_scans(["scan0001", "scan0002", "scan0003"])

    view.remove_scans(["scan0001", "scan0003"])

    assert view.scan_names() == ["scan0002"]


def test_remove_scans_ignores_names_that_are_not_listed(view):
    view.add_scans(["scan0001"])

    view.remove_scans(["scan0099"])

    assert view.scan_names() == ["scan0001"]


def test_removed_scan_can_be_added_again(view):
    view.add_scans(["scan0001"])
    view.remove_scans(["scan0001"])
    view.add_scans(["scan0001"])

    assert view.scan_names() == ["scan0001"]


def test_selected_scan_names_reports_the_selection(view):
    view.add_scans(["scan0001", "scan0002", "scan0003"])
    _select(view, "scan0001", "scan0003")

    assert view.selected_scan_names() == ["scan0001", "scan0003"]


def test_clear_scans_disables_the_remove_button(view):
    view.add_scans(["scan0001"])
    _select(view, "scan0001")

    view.clear_scans()

    assert not view.remove_scans_btn.isEnabled()


# ---------------------------------------------------------------------------
# CombineView — scan list size
# ---------------------------------------------------------------------------


def test_scan_list_is_capped_at_half_its_natural_height(view):
    assert view.scan_list.maximumHeight() == QListWidget().sizeHint().height() // 2


# ---------------------------------------------------------------------------
# CombineView — rebinning
# ---------------------------------------------------------------------------


def test_rebin_tabs_are_equal_step_tolerance_and_histogram(view):
    labels = [view.rebin_tabs.tabText(i) for i in range(view.rebin_tabs.count())]

    assert labels == ["Equal Step", "Tolerance", "Histogram"]


def test_rebin_mode_defaults_to_equal_step(view):
    assert view.rebin_mode() is RebinMode.EQUAL_STEP


def test_rebin_mode_follows_the_visible_tab(view):
    view.rebin_tabs.setCurrentIndex(1)

    assert view.rebin_mode() is RebinMode.TOLERANCE


def test_rebin_mode_is_histogram_on_the_histogram_tab(view):
    view.rebin_tabs.setCurrentIndex(2)

    assert view.rebin_mode() is RebinMode.HISTOGRAM


def _histogram_edits(view):
    return (
        view.histogram_x_start_edit,
        view.histogram_x_stop_edit,
        view.histogram_x_step_edit,
        view.histogram_y_start_edit,
        view.histogram_y_stop_edit,
        view.histogram_y_step_edit,
    )


def _histogram_values(view):
    return (
        view.histogram_x_start(),
        view.histogram_x_stop(),
        view.histogram_x_step(),
        view.histogram_y_start(),
        view.histogram_y_stop(),
        view.histogram_y_step(),
    )


def test_histogram_default_field_values(view):
    """Start and Stop default to none on both axes, so the grid's edges are read off the data."""
    assert view.histogram_z_axis() == "e"
    assert _histogram_values(view) == ("", "", "0.02", "", "", "0.02")


def test_histogram_blank_edges_show_none_as_a_placeholder(view):
    for edit in (
        view.histogram_x_start_edit,
        view.histogram_x_stop_edit,
        view.histogram_y_start_edit,
        view.histogram_y_stop_edit,
    ):
        assert edit.placeholderText() == "None"


def test_histogram_fields_report_what_was_typed(view):
    view.histogram_z_axis_combo.setEditText("qh")
    typed = ("1.5", "2.5", "0.05", "-1", "8", "0.5")
    for edit, text in zip(_histogram_edits(view), typed):
        edit.setText(text)

    assert view.histogram_z_axis() == "qh"
    assert _histogram_values(view) == typed


def test_histogram_tab_holds_only_its_own_fields(view):
    histogram_page = view.rebin_tabs.widget(2)

    for field in (view.histogram_z_axis_combo, *_histogram_edits(view)):
        assert field.parent() is histogram_page
    for field in (view.rebin_tolerance_edit, view.rebin_start_edit, view.rebin_stop_edit, view.rebin_step_edit):
        assert field.parent() is not histogram_page


def _row_texts(row):
    texts = []
    for i in range(row.count()):
        widget = row.itemAt(i).widget()
        if isinstance(widget, QLabel):
            texts.append(widget.text())
        elif widget is not None:
            texts.append(widget)
    return texts


def test_histogram_tab_stacks_z_axis_then_x_then_y(view):
    page_layout = view.rebin_tabs.widget(2).layout()
    rows = [_row_texts(page_layout.itemAt(i).layout()) for i in range(page_layout.count())]

    assert rows == [
        ["Z-Axis", view.histogram_z_axis_combo],
        [
            "X:",
            "Start",
            view.histogram_x_start_edit,
            "Stop",
            view.histogram_x_stop_edit,
            "Step",
            view.histogram_x_step_edit,
        ],
        [
            "Y:",
            "Start",
            view.histogram_y_start_edit,
            "Stop",
            view.histogram_y_stop_edit,
            "Step",
            view.histogram_y_step_edit,
        ],
    ]
    assert isinstance(page_layout, QVBoxLayout)


def test_rebin_note_explains_the_histogram(view):
    view.rebin_tabs.setCurrentIndex(2)

    assert view.rebin_note.text() == view.REBIN_NOTES[RebinMode.HISTOGRAM]


def test_equal_step_default_field_values(view):
    """Start and Stop default to none, so the grid's edges are read off the data."""
    assert view.rebin_start() == ""
    assert view.rebin_stop() == ""
    assert view.rebin_step() == "0.02"


def test_rebin_note_explains_the_equal_step_bin_edges(view):
    view.rebin_tabs.setCurrentIndex(0)

    assert view.rebin_note.text() == (
        "*When Start/Stop are set to None, bin edges are [min - step/2, max + step/2].\n"
        "*When Start/Stop are specified, they define bin edge."
    )


def test_rebin_note_explains_the_tolerance_bin_edges(view):
    view.rebin_tabs.setCurrentIndex(1)

    assert view.rebin_note.text() == "*Bin edges are [min - step/2, max + step/2]."


def test_rebin_note_follows_the_visible_tab_back_and_forth(view):
    view.rebin_tabs.setCurrentIndex(1)
    view.rebin_tabs.setCurrentIndex(0)

    assert view.rebin_note.text() == view.REBIN_NOTES[RebinMode.EQUAL_STEP]


def test_every_rebin_tab_has_a_note(view):
    assert set(view.REBIN_NOTES) == set(view._REBIN_MODE_BY_TAB)


def test_rebin_note_sits_directly_under_the_rebin_tabs(view):
    box_layout = view.rebin_tabs.parentWidget().layout()

    assert box_layout.indexOf(view.rebin_note) == box_layout.indexOf(view.rebin_tabs) + 1


def test_rebin_note_is_grey(view):
    assert "color: gray" in view.rebin_note.styleSheet()


def test_equal_step_blank_edges_show_none_as_a_placeholder(view):
    assert view.rebin_start_edit.placeholderText() == "None"
    assert view.rebin_stop_edit.placeholderText() == "None"


def test_tolerance_default_field_value(view):
    assert view.rebin_tolerance() == "0.02"


def test_equal_step_fields_report_what_was_typed(view):
    view.rebin_start_edit.setText("-1.5")
    view.rebin_stop_edit.setText("3.5")
    view.rebin_step_edit.setText("0.1")

    assert (view.rebin_start(), view.rebin_stop(), view.rebin_step()) == ("-1.5", "3.5", "0.1")


def test_tolerance_field_reports_what_was_typed(view):
    view.rebin_tolerance_edit.setText("0.005")

    assert view.rebin_tolerance() == "0.005"


def test_tolerance_tab_offers_no_grid_edges(view):
    """start/finish belong to an equal-step grid; CombineOps refuses them for a tolerance rebin."""
    tolerance_page = view.rebin_tabs.widget(1)

    assert view.rebin_tolerance_edit.parent() is tolerance_page
    assert view.rebin_start_edit.parent() is not tolerance_page
    assert view.rebin_stop_edit.parent() is not tolerance_page
    assert view.rebin_step_edit.parent() is not tolerance_page


# ---------------------------------------------------------------------------
# CombineView — axis, normalization and multiplication
# ---------------------------------------------------------------------------


def test_axis_default_values(view):
    assert view.y_axis() == "detector"
    assert view.x_axis() == "s2"


def test_normalization_defaults_to_mcu(view):
    assert view.normalization() == "mcu"


def test_multiplication_defaults_to_one(view):
    assert view.multiplication() == "1"


def test_column_fields_accept_typed_names(view):
    view.y_axis_combo.setEditText("det_2")
    view.x_axis_combo.setEditText("en")
    view.normalization_combo.setEditText("monitor")

    assert (view.y_axis(), view.x_axis(), view.normalization()) == ("det_2", "en", "monitor")


def test_multiplication_reports_what_was_typed(view):
    view.multiplication_edit.setText("2.5")

    assert view.multiplication() == "2.5"


def test_set_axes_shows_the_named_columns(view):
    view.set_axes("en", "det_2")

    assert (view.x_axis(), view.y_axis()) == ("en", "det_2")


def test_set_axes_keeps_the_index_in_step_with_the_text(view):
    view.set_column_options(["qh", "detector"])

    view.set_axes("qh", "detector")

    assert view.x_axis_combo.currentIndex() == 0
    assert view.y_axis_combo.currentIndex() == 1


def test_set_axes_leaves_the_normalization_alone(view):
    view.set_axes("qh", "detector")

    assert view.normalization() == "mcu"
    assert view.histogram_z_axis() == "e"


def test_set_column_options_fills_every_column_dropdown(view):
    view.set_column_options(["s2", "detector", "monitor"])

    for combo in (view.y_axis_combo, view.x_axis_combo, view.normalization_combo, view.histogram_z_axis_combo):
        items = [combo.itemText(i) for i in range(combo.count())]
        assert items == ["s2", "detector", "monitor"]


def test_set_column_options_keeps_the_names_already_shown(view):
    view.normalization_combo.setEditText("monitor")

    view.set_column_options(["s2", "detector", "monitor", "mcu"])

    assert view.y_axis() == "detector"
    assert view.x_axis() == "s2"
    assert view.normalization() == "monitor"


def test_set_column_options_keeps_a_name_the_scan_does_not_have(view):
    """A column can be named before the scan carrying it is focused, so a relist must not drop it."""
    view.y_axis_combo.setEditText("det_hi_res")

    view.set_column_options(["s2", "detector"])

    assert view.y_axis() == "det_hi_res"


def test_set_column_options_replaces_a_previous_list(view):
    view.set_column_options(["s2", "detector"])
    view.set_column_options(["en", "monitor"])

    items = [view.x_axis_combo.itemText(i) for i in range(view.x_axis_combo.count())]
    assert items == ["en", "monitor"]


# ---------------------------------------------------------------------------
# CombineView — Append / Combine
# ---------------------------------------------------------------------------


def test_run_button_labels(view):
    assert view.append_btn.text() == "Append"
    assert view.combine_btn.text() == "Combine"


def test_run_buttons_disabled_while_no_scans_are_listed(view):
    assert not view.append_btn.isEnabled()
    assert not view.combine_btn.isEnabled()


def test_run_buttons_enabled_once_a_scan_is_listed(view):
    view.add_scans(["scan0001"])

    assert view.append_btn.isEnabled()
    assert view.combine_btn.isEnabled()


def test_run_buttons_disabled_again_when_the_last_scan_is_removed(view):
    view.add_scans(["scan0001"])
    view.remove_scans(["scan0001"])

    assert not view.append_btn.isEnabled()
    assert not view.combine_btn.isEnabled()


def test_run_buttons_disabled_after_clear_scans(view):
    view.add_scans(["scan0001", "scan0002"])
    view.clear_scans()

    assert not view.append_btn.isEnabled()
    assert not view.combine_btn.isEnabled()


def test_append_button_emits_append_clicked(view, qtbot):
    view.add_scans(["scan0001"])

    with qtbot.waitSignal(view.append_clicked, timeout=1000):
        view.append_btn.click()


def test_combine_button_emits_combine_clicked(view, qtbot):
    view.add_scans(["scan0001"])

    with qtbot.waitSignal(view.combine_clicked, timeout=1000):
        view.combine_btn.click()


def _layout_index_of_widget(view, widget):
    """Position of the top-level row holding widget, whether it sits in the layout or a sub-layout."""
    layout = view.layout()
    for index in range(layout.count()):
        item = layout.itemAt(index)
        if item.widget() is widget:
            return index
        sub = item.layout()
        if sub is not None and any(sub.itemAt(i).widget() is widget for i in range(sub.count())):
            return index
    raise AssertionError(f"{widget} is not in the top-level layout")


def test_run_buttons_sit_directly_below_the_rebinning_box(view):
    box = next(child for child in view.findChildren(QGroupBox) if child.title() == "Rebinning")

    assert _layout_index_of_widget(view, view.append_btn) == _layout_index_of_widget(view, box) + 1


def test_append_and_combine_share_one_row(view):
    assert _layout_index_of_widget(view, view.append_btn) == _layout_index_of_widget(view, view.combine_btn)


def test_rebin_tab_fields_share_one_row(view):
    """A tab one line tall, rather than a form row per field - except Histogram, which stacks its axes."""
    for index in range(view.rebin_tabs.count()):
        if view.rebin_tabs.tabText(index) == "Histogram":
            continue
        assert isinstance(view.rebin_tabs.widget(index).layout(), QHBoxLayout)
