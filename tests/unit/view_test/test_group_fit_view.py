"""Tests for GroupFitWindow."""

import pytest
from qtpy.QtCore import QItemSelectionModel
from qtpy.QtWidgets import QAbstractItemView, QLabel, QSplitter

from tavi.frontend.view.fitting_view import FittingView
from tavi.frontend.view.group_fit_view import NO_SCAN_TEXT, GroupFitWindow
from tavi.library.data.fit_entry import FitCurve, FitData, FitResultSummary, FitSpec, ParamField, PeakField, PeakResult
from tavi.library.data.group_fit_entry import GroupFitOutcome
from tavi.library.data.scan import UUID


def make_param(value="", fixed=False, minimum="", maximum="") -> ParamField:
    return ParamField(value=value, fixed=fixed, minimum=minimum, maximum=maximum)


def make_spec(range_min="0", range_max="10", center="0", **overrides) -> FitSpec:
    fields = dict(
        range_min=range_min,
        range_max=range_max,
        background="Linear",
        background_constant=make_param("0"),
        background_slope=make_param("0"),
        peaks=[PeakField(shape="Gaussian", amplitude=make_param("1"), center=make_param(center), fwhm=make_param("1"))],
    )
    fields.update(overrides)
    return FitSpec(**fields)


def make_scans(*names: str) -> list[tuple[UUID, str]]:
    return [(UUID(value=f"uuid-{name}"), name) for name in names]


def list_scans(window, scans, spec=None) -> None:
    """
    List ``scans``; in an empty window each starts from ``spec``, typed into the panel first as a user would.

    A window that already lists scans is left as it is, so new scans start from its first scan's spec.
    """
    if window.scan_list.count() == 0:
        window.fitting_view.set_fit_spec(spec if spec is not None else make_spec())
    window.set_scans(scans)


@pytest.fixture
def window(qtbot):
    w = GroupFitWindow(FittingView())
    qtbot.addWidget(w)
    return w


def select(window: GroupFitWindow, row: int) -> None:
    """A plain click: pick ``row`` alone and make it current."""
    window.scan_list.setCurrentRow(row, QItemSelectionModel.ClearAndSelect)


# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------


def test_window_starts_with_no_scans(window):
    assert window.scan_list.count() == 0
    assert window.current_scan_uuid() is None


def test_panel_offers_no_undo_or_redo_fit(window):
    window.show()

    assert not window.fitting_view.undo_fit_btn.isVisible()
    assert not window.fitting_view.redo_fit_btn.isVisible()


def test_panel_hides_the_1d_tabs_add_to_group_fit_button(window):
    window.show()

    assert not window.fitting_view.group_fit_btn.isVisible()


def test_panel_fit_button_is_group_fit(window):
    window.show()

    assert window.fitting_view.perform_fit_btn.text() == "Group Fit"
    assert window.fitting_view.perform_fit_btn.isVisible()


def test_group_fit_button_emits_group_fit(window):
    window.set_scans(make_scans("s1"))
    received = []
    window.hookup_group_fit_signal(lambda: received.append(True))

    window.fitting_view.perform_fit_btn.click()

    assert received == [True]


def test_panel_is_disabled_until_a_scan_is_picked(window):
    assert not window.fitting_view.isEnabled()


# ---------------------------------------------------------------------------
# set_scans
# ---------------------------------------------------------------------------


def test_set_scans_lists_names_in_order(window):
    window.set_scans(make_scans("s1", "s2", "s3"))

    assert [window.scan_list.item(row).text() for row in range(3)] == ["s1", "s2", "s3"]


def test_set_scans_picks_the_first_scan(window):
    scans = make_scans("s1", "s2")

    window.set_scans(scans)

    assert window.current_scan_uuid() == scans[0][0]
    assert window.fitting_view.isEnabled()


def test_set_scans_into_an_empty_window_seeds_from_the_panel(window):
    window.fitting_view.min_edit.setText("-7")

    window.set_scans(make_scans("s1", "s2"))

    assert [spec.range_min for spec in window.fit_specs().values()] == ["-7", "-7"]


def test_set_scans_keeps_specs_of_scans_still_listed(window):
    s1, s2, s3 = make_scans("s1", "s2", "s3")
    list_scans(window, [s1, s2])
    select(window, 1)
    window.fitting_view.min_edit.setText("-4")

    window.set_scans([s1, s2, s3])

    specs = window.fit_specs()
    assert [specs[uuid].range_min for uuid, _ in (s1, s2, s3)] == ["0", "-4", "0"]


def test_set_scans_seeds_a_new_scan_from_the_first_scans_spec(window):
    s1, s2, s3 = make_scans("s1", "s2", "s3")
    list_scans(window, [s1, s2], make_spec(range_min="-6"))
    select(window, 1)
    # The panel shows s2, so a seed taken from the panel would carry this edit; the group model doesn't.
    window.fitting_view.max_edit.setText("4")

    window.set_scans([s1, s2, s3])

    spec = window.fit_specs()[s3[0]]
    assert (spec.range_min, spec.range_max) == ("-6", "10")


def test_set_scans_keeps_the_current_scan_on_screen(window):
    s1, s2, s3 = make_scans("s1", "s2", "s3")
    list_scans(window, [s1, s2])
    select(window, 1)

    window.set_scans([s3, s2, s1])

    assert window.current_scan_uuid() == s2[0]
    assert window.scan_list.currentRow() == 1


def test_set_scans_drops_scans_no_longer_listed(window):
    s1, s2 = make_scans("s1", "s2")
    list_scans(window, [s1, s2])
    select(window, 1)

    window.set_scans([s1])

    assert list(window.fit_specs()) == [s1[0]]
    assert window.current_scan_uuid() == s1[0]


def test_set_scans_empty_disables_the_panel(window):
    window.set_scans(make_scans("s1"))

    window.set_scans([])

    assert window.current_scan_uuid() is None
    assert not window.fitting_view.isEnabled()
    assert window.fit_specs() == {}


# ---------------------------------------------------------------------------
# Switching scans
# ---------------------------------------------------------------------------


def test_switching_scans_shows_each_scans_own_spec(window):
    list_scans(window, make_scans("s1", "s2"))
    window.fitting_view.min_edit.setText("-1")
    select(window, 1)
    window.fitting_view.min_edit.setText("-2")

    select(window, 0)
    assert window.fitting_view.min_edit.text() == "-1"
    select(window, 1)
    assert window.fitting_view.min_edit.text() == "-2"


def test_switching_scans_restores_peak_count_and_shape(window):
    list_scans(window, make_scans("s1", "s2"))
    select(window, 1)
    window.fitting_view.num_peaks_spin.setValue(2)
    window.fitting_view.peak_panels[1].peak_shape_combo.setCurrentText("Cosine")

    select(window, 0)
    assert len(window.fitting_view.peak_panels) == 1

    select(window, 1)
    assert [panel.peak_shape_combo.currentText() for panel in window.fitting_view.peak_panels] == [
        "Gaussian",
        "Cosine",
    ]


def test_fit_specs_include_unstored_edits_on_screen(window):
    scans = make_scans("s1")
    list_scans(window, scans)

    window.fitting_view.max_edit.setText("42")

    assert window.fit_specs()[scans[0][0]].range_max == "42"


def test_picking_a_scan_emits_its_uuid(window):
    scans = make_scans("s1", "s2")
    window.set_scans(scans)
    received = []
    window.hookup_scan_selected_signal(received.append)

    select(window, 1)

    assert received == [scans[1][0]]


# ---------------------------------------------------------------------------
# Plot placeholder
# ---------------------------------------------------------------------------


def test_window_has_a_canvas_with_one_empty_axes(window):
    assert window.canvas.figure is window.figure
    assert window.figure.axes == [window.axes]
    assert not window.axes.lines


def test_canvas_sits_between_the_scan_list_and_the_panel(window):
    splitter = window.findChild(QSplitter)

    assert splitter.count() == 3
    assert splitter.widget(0).isAncestorOf(window.scan_list)
    assert splitter.widget(1).isAncestorOf(window.canvas)
    assert splitter.widget(2).isAncestorOf(window.fitting_view)


# ---------------------------------------------------------------------------
# Remove Scan
# ---------------------------------------------------------------------------


def test_remove_scan_is_disabled_until_a_scan_is_picked(window):
    assert not window.remove_scan_btn.isEnabled()

    window.set_scans(make_scans("s1"))

    assert window.remove_scan_btn.isEnabled()


def test_remove_scan_drops_the_picked_scan_and_its_spec(window):
    s1, s2, s3 = make_scans("s1", "s2", "s3")
    list_scans(window, [s1, s2, s3])
    select(window, 1)

    window.remove_scan_btn.click()

    assert [window.scan_list.item(row).text() for row in range(window.scan_list.count())] == ["s1", "s3"]
    assert list(window.fit_specs()) == [s1[0], s3[0]]


def test_remove_scan_moves_on_to_the_next_scan_and_shows_its_spec(window):
    s1, s2, s3 = make_scans("s1", "s2", "s3")
    list_scans(window, [s1, s2, s3])
    select(window, 2)
    window.fitting_view.min_edit.setText("-3")
    select(window, 1)
    window.fitting_view.min_edit.setText("-2")

    window.remove_scan_btn.click()

    assert window.current_scan_uuid() == s3[0]
    assert window.fitting_view.min_edit.text() == "-3"


def test_removing_the_last_listed_scan_moves_back_to_the_previous_one(window):
    s1, s2 = make_scans("s1", "s2")
    list_scans(window, [s1, s2])
    select(window, 1)

    window.remove_scan_btn.click()

    assert window.current_scan_uuid() == s1[0]


def test_removing_the_only_scan_disables_the_panel_and_the_button(window):
    window.set_scans(make_scans("s1"))

    window.remove_scan_btn.click()

    assert window.scan_list.count() == 0
    assert window.current_scan_uuid() is None
    assert window.fit_specs() == {}
    assert not window.fitting_view.isEnabled()
    assert not window.remove_scan_btn.isEnabled()


def test_remove_scan_emits_the_removed_uuid(window):
    scans = make_scans("s1", "s2")
    window.set_scans(scans)
    received = []
    window.hookup_scan_removed_signal(received.append)

    window.remove_scan_btn.click()

    assert received == [scans[0][0]]


# ---------------------------------------------------------------------------
# Fitting tab
# ---------------------------------------------------------------------------


def test_fitting_panel_sits_in_a_single_tab(window):
    assert window.fitting_tabs.count() == 1
    assert window.fitting_tabs.widget(0).isAncestorOf(window.fitting_view)


def test_fitting_tab_starts_with_no_scan(window):
    assert window.fitting_tabs.tabText(0) == NO_SCAN_TEXT


def test_fitting_tab_is_named_after_the_picked_scan(window):
    window.set_scans(make_scans("s1", "s2"))
    assert window.fitting_tabs.tabText(0) == "s1"

    select(window, 1)

    assert window.fitting_tabs.tabText(0) == "s2"


def test_fitting_tab_name_resets_once_the_last_scan_is_removed(window):
    window.set_scans(make_scans("s1"))

    window.remove_scan_btn.click()

    assert window.fitting_tabs.tabText(0) == NO_SCAN_TEXT


# ---------------------------------------------------------------------------
# add_scans
# ---------------------------------------------------------------------------


def test_add_scans_appends_after_the_listed_scans(window):
    s1, s2, s3 = make_scans("s1", "s2", "s3")
    window.set_scans([s1])

    window.add_scans([s2, s3])

    assert [window.scan_list.item(row).text() for row in range(3)] == ["s1", "s2", "s3"]


def test_add_scans_skips_scans_already_listed_and_keeps_their_spec(window):
    s1, s2, s3 = make_scans("s1", "s2", "s3")
    list_scans(window, [s1, s2])
    select(window, 1)
    window.fitting_view.min_edit.setText("-5")

    window.add_scans([s2, s3])

    specs = window.fit_specs()
    assert window.scan_list.count() == 3
    # s2 keeps its own edit; s3 joins from the group model, s1's spec.
    assert [specs[uuid].range_min for uuid, _ in (s1, s2, s3)] == ["0", "-5", "0"]


def test_add_scans_keeps_the_scan_on_screen(window):
    s1, s2, s3 = make_scans("s1", "s2", "s3")
    window.set_scans([s1, s2])
    select(window, 1)

    window.add_scans([s3])

    assert window.current_scan_uuid() == s2[0]


def test_add_scans_to_an_empty_window_picks_the_first_added(window):
    scans = make_scans("s1", "s2")

    window.add_scans(scans)

    assert window.current_scan_uuid() == scans[0][0]
    assert window.fitting_view.isEnabled()


# ---------------------------------------------------------------------------
# Fit y vs x columns
# ---------------------------------------------------------------------------


def finish_editing(edit, text):
    edit.setText(text)
    edit.editingFinished.emit()


def test_columns_row_sits_above_the_canvas(window):
    plot_layout = window.canvas.parentWidget().layout()
    columns_row = plot_layout.itemAt(0).layout()

    assert columns_row.indexOf(window.y_col_edit) < columns_row.indexOf(window.x_col_edit)
    assert plot_layout.itemAt(1).layout().indexOf(window.title_edit) >= 0
    assert plot_layout.indexOf(window.canvas) == 2


def test_columns_row_reads_fit_y_vs_x(window):
    row = window.canvas.parentWidget().layout().itemAt(0).layout()
    labels = [row.itemAt(i).widget().text() for i in range(row.count()) if isinstance(row.itemAt(i).widget(), QLabel)]

    assert labels == ["Fit y:", "vs x:"]


def test_columns_start_blank(window):
    assert window.columns() == ("", "")


def test_columns_are_read_as_x_then_y_without_surrounding_spaces(window):
    window.y_col_edit.setText(" detector ")
    window.x_col_edit.setText("s1 ")

    assert window.columns() == ("s1", "detector")


def test_editing_a_column_emits_x_and_y(window):
    received = []
    window.hookup_columns_changed_signal(lambda x, y: received.append((x, y)))

    finish_editing(window.y_col_edit, "detector")
    finish_editing(window.x_col_edit, "s1")

    assert received == [("", "detector"), ("s1", "detector")]


def test_leaving_a_column_unchanged_emits_nothing(window):
    finish_editing(window.y_col_edit, "detector")
    received = []
    window.hookup_columns_changed_signal(lambda x, y: received.append((x, y)))

    window.y_col_edit.editingFinished.emit()
    finish_editing(window.x_col_edit, "")

    assert received == []


def test_set_columns_shows_them_without_emitting(window):
    received = []
    window.hookup_columns_changed_signal(lambda x, y: received.append((x, y)))

    window.set_columns("s1", "detector")
    window.y_col_edit.editingFinished.emit()

    assert (window.x_col_edit.text(), window.y_col_edit.text()) == ("s1", "detector")
    assert received == []


def test_columns_are_shared_by_every_scan(window):
    window.set_scans(make_scans("s1", "s2"))
    finish_editing(window.y_col_edit, "detector")

    select(window, 1)

    assert window.columns() == ("", "detector")


# ---------------------------------------------------------------------------
# Overplot of the selected scans
# ---------------------------------------------------------------------------


def make_outcome(uuid, offset=0.0, components=None) -> GroupFitOutcome:
    x = [0.0, 1.0, 2.0]
    y = [offset + v for v in (1.0, 3.0, 1.0)]
    return GroupFitOutcome(
        source_scan_uuid=uuid,
        spec=make_spec(),
        data=FitData(x=x, y=y, err=[0.1, 0.1, 0.1]),
        curve=FitCurve(source_scan_uuid=uuid, scan_name=uuid.value, x=x, best_fit=y, components=components or {}),
    )


def fitted_window(window, *names, components=None):
    scans = make_scans(*names)
    list_scans(window, scans)
    window.show_outcomes(None, "detector", [make_outcome(uuid, i, components) for i, (uuid, _) in enumerate(scans)])
    return scans


def add_to_selection(window, row):
    """Ctrl-click: select ``row`` too, and make it current."""
    window.scan_list.setCurrentRow(row, QItemSelectionModel.Select)


def drawn_labels(window):
    return window.axes.get_legend_handles_labels()[1]


def test_scan_list_allows_selecting_several_scans(window):
    assert window.scan_list.selectionMode() == QAbstractItemView.ExtendedSelection


def test_nothing_is_drawn_before_a_group_fit(window):
    window.set_scans(make_scans("s1", "s2"))

    assert drawn_labels(window) == []


def test_outcomes_draw_the_selected_scans_data_and_fit(window):
    fitted_window(window, "s1", "s2")

    assert drawn_labels(window) == ["s1 fit", "s1"]


def test_selecting_more_scans_overplots_them(window):
    fitted_window(window, "s1", "s2", "s3")

    add_to_selection(window, 2)

    assert sorted(drawn_labels(window)) == ["s1", "s1 fit", "s3", "s3 fit"]
    assert len(window.axes.containers) == 2


def test_picking_a_single_scan_replaces_the_overplot(window):
    fitted_window(window, "s1", "s2", "s3")
    add_to_selection(window, 2)

    select(window, 1)

    assert sorted(drawn_labels(window)) == ["s2", "s2 fit"]


def test_the_current_scan_is_the_last_one_clicked(window):
    s1, s2 = fitted_window(window, "s1", "s2")

    add_to_selection(window, 1)

    assert window.selected_scan_uuids() == [s1[0], s2[0]]
    assert window.current_scan_uuid() == s2[0]
    assert window.fitting_tabs.tabText(0) == "s2"


def test_a_scan_keeps_its_color_as_the_selection_changes(window):
    fitted_window(window, "s1", "s2")
    select(window, 1)
    alone = window.axes.lines[0].get_color()

    add_to_selection(window, 0)

    colors = {line.get_label(): line.get_color() for line in window.axes.lines}
    assert colors["s2 fit"] == alone
    assert colors["s1 fit"] != alone


def test_a_scans_data_and_fit_share_a_color(window):
    fitted_window(window, "s1")

    (container,) = window.axes.containers
    assert container.lines[0].get_color() == window.axes.lines[-1].get_color()


def test_a_scan_without_an_outcome_is_not_drawn(window):
    scans = make_scans("s1", "s2")
    list_scans(window, scans)
    window.show_outcomes(None, "detector", [make_outcome(scans[0][0])])

    select(window, 1)

    assert drawn_labels(window) == []


def test_axes_are_labelled_with_the_fitted_columns(window):
    scans = fitted_window(window, "s1")
    assert (window.axes.get_xlabel(), window.axes.get_ylabel()) == ("default x axis", "detector")

    window.show_outcomes("s2", "detector", [make_outcome(scans[0][0])])

    assert window.axes.get_xlabel() == "s2"


def test_plot_separately_draws_each_component(window):
    fitted_window(window, "s1", components={"peak1_": [0.0, 1.0, 0.0], "bg_": [1.0, 1.0, 1.0]})
    assert "s1 peak1" not in drawn_labels(window)

    window.fitting_view.plot_sep_check.setChecked(True)

    assert {"s1 peak1", "s1 bg"} <= set(drawn_labels(window))


def test_outcomes_for_scans_no_longer_listed_are_dropped(window):
    s1, s2 = make_scans("s1", "s2")
    list_scans(window, [s1])

    window.show_outcomes(None, "detector", [make_outcome(s1[0]), make_outcome(s2[0])])
    window.add_scans([s2])
    select(window, 1)

    assert drawn_labels(window) == []


def test_show_outcomes_signal_redraws(window, qtbot):
    (s1,) = make_scans("s1")
    list_scans(window, [s1])

    with qtbot.waitSignal(window.show_outcomes_signal, timeout=1000):
        window.show_outcomes_signal.emit(None, "detector", [make_outcome(s1[0])])

    assert drawn_labels(window) == ["s1 fit", "s1"]


# ---------------------------------------------------------------------------
# Removing several scans
# ---------------------------------------------------------------------------


def test_remove_scan_drops_every_selected_scan(window):
    s1, s2, s3 = fitted_window(window, "s1", "s2", "s3")
    add_to_selection(window, 2)
    received = []
    window.hookup_scan_removed_signal(received.append)

    window.remove_scan_btn.click()

    assert [window.scan_list.item(row).text() for row in range(window.scan_list.count())] == ["s2"]
    assert received == [s1[0], s3[0]]
    assert window.current_scan_uuid() == s2[0]
    assert drawn_labels(window) == ["s2 fit", "s2"]


def test_remove_scans_keeps_the_current_scan_and_its_edits(window):
    s1, s2, s3 = make_scans("s1", "s2", "s3")
    list_scans(window, [s1, s2, s3])
    select(window, 1)
    window.fitting_view.min_edit.setText("-6")

    window.remove_scans([s1[0]])

    assert window.current_scan_uuid() == s2[0]
    assert window.fitting_view.min_edit.text() == "-6"
    assert window.fit_specs()[s2[0]].range_min == "-6"


def test_remove_scans_ignores_scans_not_listed(window):
    window.set_scans(make_scans("s1"))
    received = []
    window.hookup_scan_removed_signal(received.append)

    window.remove_scans([UUID(value="elsewhere")])

    assert window.scan_list.count() == 1
    assert received == []


def test_removing_the_current_scan_loads_the_next_scans_spec(window):
    s1, s2 = make_scans("s1", "s2")
    list_scans(window, [s1, s2])
    select(window, 1)
    window.fitting_view.min_edit.setText("-8")
    select(window, 0)

    window.remove_scans([s1[0]])

    assert window.current_scan_uuid() == s2[0]
    assert window.fitting_view.min_edit.text() == "-8"


def test_add_scans_signal_lists_the_scans_and_opens_the_window(window, qtbot):
    with qtbot.waitSignal(window.add_scans_signal, timeout=1000):
        window.add_scans_signal.emit(make_scans("s1"))

    assert window.isVisible()
    assert window.scan_list.item(0).text() == "s1"


# ---------------------------------------------------------------------------
# The first scan is the group model
# ---------------------------------------------------------------------------


def test_window_starts_with_a_blank_range_meaning_the_whole_scan(window):
    assert (window.fitting_view.min_edit.text(), window.fitting_view.max_edit.text()) == ("", "")
    assert window.fitting_view.min_edit.placeholderText() == "whole scan"


def test_editing_the_first_scan_reaches_the_others(window):
    s1, s2, s3 = make_scans("s1", "s2", "s3")
    list_scans(window, [s1, s2, s3])

    window.fitting_view.peak_panels[0].peak_table.rows[1].value_edit.setText("-0.6")
    window.fitting_view.min_edit.setText("-2")

    specs = window.fit_specs()
    assert [specs[uuid].peaks[0].center.value for uuid, _ in (s1, s2, s3)] == ["-0.6"] * 3
    assert [specs[uuid].range_min for uuid, _ in (s1, s2, s3)] == ["-2"] * 3


def test_a_field_edited_on_its_own_scan_keeps_its_value(window):
    s1, s2 = make_scans("s1", "s2")
    list_scans(window, [s1, s2])
    select(window, 1)
    window.fitting_view.peak_panels[0].peak_table.rows[1].value_edit.setText("0.4")
    select(window, 0)

    window.fitting_view.peak_panels[0].peak_table.rows[1].value_edit.setText("-0.6")
    window.fitting_view.max_edit.setText("3")

    spec = window.fit_specs()[s2[0]]
    assert spec.peaks[0].center.value == "0.4"
    assert spec.range_max == "3"


def test_a_first_scan_edit_reaches_a_fixed_flag_too(window):
    s1, s2 = make_scans("s1", "s2")
    list_scans(window, [s1, s2])

    window.fitting_view.peak_panels[0].peak_table.rows[2].fix_check.setChecked(True)

    assert window.fit_specs()[s2[0]].peaks[0].fwhm.fixed is True


def test_editing_another_scan_does_not_reach_the_rest(window):
    s1, s2, s3 = make_scans("s1", "s2", "s3")
    list_scans(window, [s1, s2, s3])
    select(window, 1)

    window.fitting_view.min_edit.setText("-9")

    specs = window.fit_specs()
    assert [specs[uuid].range_min for uuid, _ in (s1, s2, s3)] == ["0", "-9", "0"]


def test_a_first_scan_peak_count_change_reaches_scans_sharing_its_peaks(window):
    s1, s2, s3 = make_scans("s1", "s2", "s3")
    list_scans(window, [s1, s2, s3])
    select(window, 2)
    window.fitting_view.peak_panels[0].peak_shape_combo.setCurrentText("Lorentzian")
    select(window, 0)

    window.fitting_view.num_peaks_spin.setValue(2)

    specs = window.fit_specs()
    assert [peak.shape for peak in specs[s2[0]].peaks] == ["Gaussian", "Gaussian"]
    # s3 already had peaks of its own, so it keeps them.
    assert [peak.shape for peak in specs[s3[0]].peaks] == ["Lorentzian"]


def test_a_first_scan_cosine_edit_reaches_the_others(window):
    s1, s2 = make_scans("s1", "s2")
    list_scans(window, [s1, s2])
    window.fitting_view.peak_panels[0].peak_shape_combo.setCurrentText("Cosine")
    window.fit_specs()

    window.fitting_view.peak_panels[0].peak_table.rows[1].value_edit.setText("2.5")

    assert window.fit_specs()[s2[0]].peaks[0].cosine.frequency.value == "2.5"


def test_scans_added_later_start_from_the_first_scan_not_the_one_on_screen(window):
    s1, s2, s3 = make_scans("s1", "s2", "s3")
    list_scans(window, [s1, s2])
    select(window, 1)
    window.fitting_view.min_edit.setText("-9")

    window.add_scans([s3])

    assert window.fit_specs()[s3[0]].range_min == "0"


# ---------------------------------------------------------------------------
# Suggest Params.
# ---------------------------------------------------------------------------


def test_no_guess_request_without_a_scan_on_screen(window):
    assert window.get_suggest_peak_request() is None
    assert window.get_suggest_background_request() is None


def test_peak_guess_request_names_the_scan_columns_range_and_asking_peak(window):
    s1, s2 = make_scans("s1", "s2")
    list_scans(window, [s1, s2], make_spec(range_min="-2", range_max=""))
    window.set_columns("", "detector")
    select(window, 1)
    window.fitting_view.num_peaks_spin.setValue(2)
    window.fitting_view.peak_panels[1].peak_shape_combo.setCurrentText("Lorentzian")
    window.fitting_view.peak_panels[1].suggest_btn.click()

    request = window.get_suggest_peak_request()

    assert request.source_scan_uuid == s2[0]
    assert (request.x_col, request.y_col) == (None, "detector")
    assert (request.range_min, request.range_max) == ("-2", "")
    assert request.shape == "Lorentzian"


def test_background_guess_request_names_the_background(window):
    list_scans(window, make_scans("s1"))
    window.set_columns("s1", "detector")

    request = window.get_suggest_background_request()

    assert (request.background, request.x_col) == ("Linear", "s1")


def test_suggest_buttons_reach_their_hookups(window):
    list_scans(window, make_scans("s1"))
    peaks, backgrounds = [], []
    window.hookup_suggest_params_signal(lambda: peaks.append(True))
    window.hookup_suggest_background_signal(lambda: backgrounds.append(True))

    window.fitting_view.peak_panels[0].suggest_btn.click()
    window.fitting_view.suggest_background_btn.click()

    assert (peaks, backgrounds) == ([True], [True])


def test_a_peak_guess_fills_the_scan_on_screen(window):
    (s1,) = make_scans("s1")
    list_scans(window, [s1])

    window.peak_params_suggested_signal.emit(s1[0], 12.0, 0.5, 0.3)

    rows = window.fitting_view.peak_panels[0].peak_table.rows
    assert [row.value_edit.text() for row in rows] == ["12", "0.5", "0.3"]


def test_a_guess_for_a_scan_no_longer_on_screen_is_dropped(window):
    s1, s2 = make_scans("s1", "s2")
    list_scans(window, [s1, s2], make_spec(center="7"))
    select(window, 1)

    window.peak_params_suggested_signal.emit(s1[0], 12.0, 0.5, 0.3)
    window.background_params_suggested_signal.emit(s1[0], 2.0, 3.0)

    assert window.fitting_view.peak_panels[0].peak_table.rows[1].value_edit.text() == "7"
    assert window.fitting_view.background_table.rows[0].value_edit.text() == "0"


def test_a_cosine_guess_fills_the_cosine(window):
    (s1,) = make_scans("s1")
    list_scans(window, [s1])
    window.fitting_view.peak_panels[0].peak_shape_combo.setCurrentText("Cosine")

    window.cosine_params_suggested_signal.emit(s1[0], 1.0, 2.0, 3.0, 4.0)

    assert [row.value_edit.text() for row in window.fitting_view.peak_panels[0].peak_table.rows] == ["1", "2", "3", "4"]


def test_a_background_guess_fills_slope_and_intercept(window):
    (s1,) = make_scans("s1")
    list_scans(window, [s1])

    window.background_params_suggested_signal.emit(s1[0], 0.25, 4.0)

    slope_row, intercept_row = window.fitting_view.background_table.rows
    assert (slope_row.value_edit.text(), intercept_row.value_edit.text()) == ("0.25", "4")


# ---------------------------------------------------------------------------
# Summarize fit results
# ---------------------------------------------------------------------------


def make_result(center=0.0) -> FitResultSummary:
    peak = PeakResult(amplitude=1.0, amplitude_err=0.1, center=center, center_err=0.01, fwhm=1.0, fwhm_err=0.1)
    return FitResultSummary(reduced_chi_squared=1.0, peaks=[peak], background_constant=0.0, background_slope=0.0)


def summarized_window(window, *names, unfitted=()):
    """List ``names`` and show a group fit result for each, except those in ``unfitted``."""
    scans = make_scans(*names)
    list_scans(window, scans)
    outcomes = [
        make_outcome(uuid).model_copy(update={"result": None if name in unfitted else make_result(float(i))})
        for i, (uuid, name) in enumerate(scans)
    ]
    window.show_outcomes(None, "detector", outcomes)
    return scans


def summary_params(window) -> list[str]:
    combo = window.summary_param_combo
    return [combo.itemText(i) for i in range(combo.count())]


def test_summary_row_sits_just_above_the_group_fit_button(window):
    rows = window.fitting_view._bottom_rows
    row = window.summary_param_combo.parentWidget()

    assert rows.indexOf(row) == 0
    assert rows.itemAt(1).layout().indexOf(window.fitting_view.perform_fit_btn) >= 0
    labels = [label.text() for label in row.findChildren(QLabel)]
    assert labels == ["Summarize fit results:", "vs"]
    assert window.generate_summary_btn.text() == "Generate"


def test_summary_offers_nothing_before_a_group_fit(window):
    list_scans(window, make_scans("s1"))

    assert summary_params(window) == []
    assert not window.generate_summary_btn.isEnabled()
    assert window.get_summarize_request() is None


def test_summary_offers_the_fitted_parameters(window):
    summarized_window(window, "s1", "s2")

    assert summary_params(window) == list(make_result().parameters())
    assert window.generate_summary_btn.isEnabled()


def test_summary_offers_nothing_when_no_scan_was_fit(window):
    summarized_window(window, "s1", unfitted=("s1",))

    assert summary_params(window) == []
    assert not window.generate_summary_btn.isEnabled()


def test_summary_keeps_its_choice_across_refits(window):
    summarized_window(window, "s1", "s2")
    window.summary_param_combo.setCurrentText("peak1 FWHM")

    summarized_window(window, "s1", "s2")

    assert window.summary_param_combo.currentText() == "peak1 FWHM"


def test_summary_request_carries_each_fitted_scans_result_in_list_order(window):
    scans = summarized_window(window, "s1", "s2", "s3", unfitted=("s2",))
    window.summary_param_combo.setCurrentText("peak1 center")
    window.summary_x_edit.setText(" temp ")

    request = window.get_summarize_request()

    assert request.parameter == "peak1 center"
    assert request.x_col == "temp"
    assert list(request.results) == [scans[0][0], scans[2][0]]
    assert request.results[scans[2][0]].peaks[0].center == 2.0


def test_summary_empties_once_every_scan_is_removed(window):
    scans = summarized_window(window, "s1")

    window.remove_scans([scans[0][0]])

    assert summary_params(window) == []
    assert not window.generate_summary_btn.isEnabled()


# ---------------------------------------------------------------------------
# Each scan's fit result in the panel
# ---------------------------------------------------------------------------


def peak_rows(window):
    return window.fitting_view.peak_panels[0].peak_table.rows


def panel_center(window) -> tuple[str, str]:
    row = peak_rows(window)[1]
    return row.value_edit.text(), row.std_edit.text()


def test_picking_a_fitted_scan_shows_its_own_result(window):
    summarized_window(window, "s1", "s2")

    select(window, 1)
    assert panel_center(window) == ("1", "0.01")
    assert window.fitting_view.chi2_edit.text() == "1"

    select(window, 0)
    assert panel_center(window) == ("0", "0.01")


def test_outcomes_arriving_show_the_current_scans_result(window):
    list_scans(window, make_scans("s1", "s2"))
    select(window, 1)

    summarized_window(window, "s1", "s2")

    assert panel_center(window) == ("1", "0.01")


def test_a_scan_left_out_of_the_fit_shows_its_spec_alone(window):
    summarized_window(window, "s1", "s2", unfitted=("s2",))

    select(window, 1)

    assert panel_center(window) == ("0", "")
    assert window.fitting_view.chi2_edit.text() == ""


def test_fitted_values_on_display_never_become_starting_values(window):
    scans = summarized_window(window, "s1", "s2")
    select(window, 1)
    select(window, 0)

    specs = window.fit_specs()

    assert [specs[uuid].peaks[0].center.value for uuid, _ in scans] == ["0", "0"]
    assert [specs[uuid].peaks[0].amplitude.value for uuid, _ in scans] == ["1", "1"]


def test_an_edit_over_a_shown_result_is_stored_alone(window):
    scans = summarized_window(window, "s1", "s2")
    select(window, 1)
    peak_rows(window)[2].value_edit.setText("0.7")

    select(window, 0)
    spec = window.fit_specs()[scans[1][0]]

    assert spec.peaks[0].fwhm.value == "0.7"
    # The center on display was the fit's, not the user's - the stored start stays.
    assert spec.peaks[0].center.value == "0"


def test_an_edit_made_while_the_fit_ran_survives_its_outcomes(window):
    scans = make_scans("s1", "s2")
    list_scans(window, scans)
    select(window, 1)
    peak_rows(window)[2].value_edit.setText("0.7")

    summarized_window(window, "s1", "s2")

    assert window.fit_specs()[scans[1][0]].peaks[0].fwhm.value == "0.7"


# ---------------------------------------------------------------------------
# Hide Label
# ---------------------------------------------------------------------------


def test_hide_label_starts_unchecked_with_the_legend_shown(window):
    fitted_window(window, "s1")
    select(window, 0)

    assert not window.hide_label_check.isChecked()
    assert window.axes.get_legend() is not None


def test_hide_label_drops_the_legend_and_keeps_the_data(window):
    fitted_window(window, "s1")
    select(window, 0)

    window.hide_label_check.setChecked(True)

    assert window.axes.get_legend() is None
    assert len(window.axes.lines) > 0


def test_hide_label_stays_hidden_as_the_selection_changes(window):
    fitted_window(window, "s1", "s2")
    window.hide_label_check.setChecked(True)

    window.scan_list.selectAll()

    assert window.axes.get_legend() is None


def test_unchecking_hide_label_restores_the_legend(window):
    fitted_window(window, "s1")
    select(window, 0)
    window.hide_label_check.setChecked(True)

    window.hide_label_check.setChecked(False)

    assert window.axes.get_legend() is not None


# ---------------------------------------------------------------------------
# Set Title
# ---------------------------------------------------------------------------


def test_title_edit_announces_its_stripped_text_once_per_change(window):
    received = []
    window.title_changed.connect(received.append)

    for text in ("  temp ", "temp", ""):
        window.title_edit.setText(text)
        window.title_edit.editingFinished.emit()

    assert received == ["temp", ""]


def test_titles_relabel_the_legend_but_not_the_scan_list(window):
    (s1, _), (s2, _) = fitted_window(window, "s1", "s2")
    window.scan_list.selectAll()

    window.set_titles([(s1, "1.5 K")])

    assert sorted(drawn_labels(window)) == ["1.5 K", "1.5 K fit", "s2", "s2 fit"]
    assert [window.scan_list.item(row).text() for row in range(2)] == ["s1", "s2"]


def test_titles_relabel_components_too(window):
    (s1, _), _ = fitted_window(window, "s1", "s2", components={"p1_": [1.0, 2.0, 1.0]})
    window.fitting_view.plot_sep_check.setChecked(True)
    select(window, 0)

    window.set_titles([(s1, "MnO")])

    assert sorted(drawn_labels(window)) == ["MnO", "MnO fit", "MnO p1"]


def test_empty_titles_restore_the_scan_names(window):
    (s1, _), _ = fitted_window(window, "s1", "s2")
    select(window, 0)
    window.set_titles([(s1, "MnO")])

    window.set_titles([])

    assert sorted(drawn_labels(window)) == ["s1", "s1 fit"]


def test_a_removed_scan_takes_its_title_with_it(window):
    (s1, _), (s2, _) = fitted_window(window, "s1", "s2")
    window.set_titles([(s1, "MnO"), (s2, "BiFeO3")])

    window.remove_scans([s1])

    assert window._titles == {s2: "BiFeO3"}


def test_titles_for_scans_not_listed_are_ignored(window):
    (s1, _), _ = fitted_window(window, "s1", "s2")

    window.set_titles([(s1, "MnO"), (UUID(value="gone"), "x")])

    assert window._titles == {s1: "MnO"}
