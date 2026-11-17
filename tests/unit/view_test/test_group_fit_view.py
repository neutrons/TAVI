"""Tests for GroupFitWindow."""

import pytest
from qtpy.QtWidgets import QSplitter

from tavi.frontend.view.fitting_view import FittingView
from tavi.frontend.view.group_fit_view import NO_SCAN_TEXT, GroupFitWindow
from tavi.library.data.fit_entry import FitSpec, ParamField, PeakField
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


@pytest.fixture
def window(qtbot):
    w = GroupFitWindow(FittingView())
    qtbot.addWidget(w)
    return w


def select(window: GroupFitWindow, row: int) -> None:
    window.scan_list.setCurrentRow(row)


# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------


def test_window_starts_with_no_scans(window):
    assert window.scan_list.count() == 0
    assert window.current_scan_uuid() is None


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


def test_set_scans_seeds_every_scan_from_the_template(window):
    scans = make_scans("s1", "s2")

    window.set_scans(scans, template=make_spec(range_min="-3"))

    assert [spec.range_min for spec in window.fit_specs().values()] == ["-3", "-3"]
    assert window.fitting_view.min_edit.text() == "-3"


def test_set_scans_without_template_seeds_from_the_panel(window):
    window.fitting_view.min_edit.setText("-7")

    window.set_scans(make_scans("s1", "s2"))

    assert [spec.range_min for spec in window.fit_specs().values()] == ["-7", "-7"]


def test_set_scans_keeps_specs_of_scans_still_listed(window):
    s1, s2, s3 = make_scans("s1", "s2", "s3")
    window.set_scans([s1, s2], template=make_spec())
    select(window, 1)
    window.fitting_view.min_edit.setText("-4")

    window.set_scans([s1, s2, s3], template=make_spec(range_min="9"))

    specs = window.fit_specs()
    assert [specs[uuid].range_min for uuid, _ in (s1, s2, s3)] == ["0", "-4", "9"]


def test_set_scans_keeps_the_current_scan_on_screen(window):
    s1, s2, s3 = make_scans("s1", "s2", "s3")
    window.set_scans([s1, s2], template=make_spec())
    select(window, 1)

    window.set_scans([s3, s2, s1])

    assert window.current_scan_uuid() == s2[0]
    assert window.scan_list.currentRow() == 1


def test_set_scans_drops_scans_no_longer_listed(window):
    s1, s2 = make_scans("s1", "s2")
    window.set_scans([s1, s2], template=make_spec())
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
    window.set_scans(make_scans("s1", "s2"), template=make_spec())
    window.fitting_view.min_edit.setText("-1")
    select(window, 1)
    window.fitting_view.min_edit.setText("-2")

    select(window, 0)
    assert window.fitting_view.min_edit.text() == "-1"
    select(window, 1)
    assert window.fitting_view.min_edit.text() == "-2"


def test_switching_scans_restores_peak_count_and_shape(window):
    window.set_scans(make_scans("s1", "s2"), template=make_spec())
    window.fitting_view.num_peaks_spin.setValue(2)
    window.fitting_view.peak_panels[1].peak_shape_combo.setCurrentText("Cosine")

    select(window, 1)
    assert len(window.fitting_view.peak_panels) == 1

    select(window, 0)
    assert [panel.peak_shape_combo.currentText() for panel in window.fitting_view.peak_panels] == [
        "Gaussian",
        "Cosine",
    ]


def test_fit_specs_include_unstored_edits_on_screen(window):
    scans = make_scans("s1")
    window.set_scans(scans, template=make_spec())

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
    window.set_scans([s1, s2, s3], template=make_spec())
    select(window, 1)

    window.remove_scan_btn.click()

    assert [window.scan_list.item(row).text() for row in range(window.scan_list.count())] == ["s1", "s3"]
    assert list(window.fit_specs()) == [s1[0], s3[0]]


def test_remove_scan_moves_on_to_the_next_scan_and_shows_its_spec(window):
    s1, s2, s3 = make_scans("s1", "s2", "s3")
    window.set_scans([s1, s2, s3], template=make_spec())
    select(window, 2)
    window.fitting_view.min_edit.setText("-3")
    select(window, 1)
    window.fitting_view.min_edit.setText("-2")

    window.remove_scan_btn.click()

    assert window.current_scan_uuid() == s3[0]
    assert window.fitting_view.min_edit.text() == "-3"


def test_removing_the_last_listed_scan_moves_back_to_the_previous_one(window):
    s1, s2 = make_scans("s1", "s2")
    window.set_scans([s1, s2], template=make_spec())
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
    s1, s2 = make_scans("s1", "s2")
    window.set_scans([s1], template=make_spec())
    window.fitting_view.min_edit.setText("-5")

    window.add_scans([s1, s2], template=make_spec(range_min="8"))

    specs = window.fit_specs()
    assert window.scan_list.count() == 2
    assert [specs[uuid].range_min for uuid, _ in (s1, s2)] == ["-5", "8"]


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
