"""Tests for MainWindow and TaviView."""

from unittest.mock import patch

import pytest
from qtpy.QtWidgets import QMenuBar, QMessageBox, QSplitter, QTabWidget

from tavi import __version__
from tavi.frontend.view.data_file_view import DataFileView
from tavi.frontend.view.filter_view import FilterView
from tavi.frontend.view.main_view import MainWindow, TaviView
from tavi.frontend.view.plotter_view import Plot1DView
from tavi.frontend.view.project_view import ProjectView


@pytest.fixture()
def built_view(qtbot):
    """TaviView with build_ui() called using real sub-views."""
    view = TaviView()
    qtbot.addWidget(view)
    view.build_ui(
        project_view=ProjectView(),
        plot_view=Plot1DView(),
        data_file_view=DataFileView(),
        filter_view=FilterView(),
    )
    return view


# ---------------------------------------------------------------------------
# MainWindow
# ---------------------------------------------------------------------------


def test_main_window_instantiates(qtbot):
    w = MainWindow()
    qtbot.addWidget(w)


def test_main_window_has_project_load_view(qtbot):
    w = MainWindow()
    qtbot.addWidget(w)
    assert isinstance(w.load_view, ProjectView)


def test_main_window_load_view_parent_is_window(qtbot):
    w = MainWindow()
    qtbot.addWidget(w)
    assert w.load_view.parent() is w


def test_handle_help_calls_help_function(qtbot):
    w = MainWindow()
    qtbot.addWidget(w)

    with patch("tavi.frontend.view.main_view.help_function") as mock_help:
        w.handle_help()
        mock_help.assert_called_once_with(context="tavi_View")


# ---------------------------------------------------------------------------
# TaviView — construction
# ---------------------------------------------------------------------------


def test_tavi_view_instantiates(qtbot):
    view = TaviView()
    qtbot.addWidget(view)


def test_tavi_view_window_title_contains_version(qtbot):
    view = TaviView()
    qtbot.addWidget(view)
    assert __version__ in view.windowTitle()


def test_tavi_view_force_closing_false_by_default(qtbot):
    view = TaviView()
    qtbot.addWidget(view)
    assert view._force_closing is False


def test_tavi_view_central_widget_is_splitter(built_view):
    assert isinstance(built_view.centralWidget(), QSplitter)


def test_tavi_view_splitter_has_three_panels(built_view):
    assert built_view.centralWidget().count() == 3


# ---------------------------------------------------------------------------
# TaviView — child widget accessors
# ---------------------------------------------------------------------------


def test_tavi_view_has_project_view(built_view):
    assert isinstance(built_view.project_view, ProjectView)


def test_tavi_view_has_plotter_view(built_view):
    assert isinstance(built_view.plotter_view, Plot1DView)


def test_data_file_view_title_changed_updates_tab_text(built_view):
    tabs = built_view.data_file_view.parent()
    while tabs is not None and not isinstance(tabs, QTabWidget):
        tabs = tabs.parent()

    built_view.data_file_view.set_title("Data File (my_scan)")

    assert tabs.tabText(0) == "Data File (my_scan)"


# ---------------------------------------------------------------------------
# TaviView — install_menu_bar
# ---------------------------------------------------------------------------


def test_install_menu_bar_sets_menu_bar(qtbot):
    view = TaviView()
    qtbot.addWidget(view)

    menu_bar = QMenuBar()
    view.install_menu_bar(menu_bar)

    assert view.menuBar() is menu_bar


def test_install_menu_bar_disables_native_mode(qtbot):
    view = TaviView()
    qtbot.addWidget(view)

    menu_bar = QMenuBar()
    view.install_menu_bar(menu_bar)

    assert not menu_bar.isNativeMenuBar()


# ---------------------------------------------------------------------------
# TaviView — closeEvent
# ---------------------------------------------------------------------------


def test_close_event_emits_exit_requested_when_not_force_closing(qtbot):
    view = TaviView()
    qtbot.addWidget(view)
    view.show()

    with qtbot.waitSignal(view.exit_requested, timeout=1000):
        view.close()


def test_close_event_keeps_window_visible_when_not_force_closing(qtbot):
    view = TaviView()
    qtbot.addWidget(view)
    view.show()

    view.close()

    assert view.isVisible()


def test_close_event_accepted_when_force_closing(qtbot):
    view = TaviView()
    qtbot.addWidget(view)
    view.show()

    view._force_closing = True
    view.close()

    assert not view.isVisible()


# ---------------------------------------------------------------------------
# TaviView — force_close
# ---------------------------------------------------------------------------


def test_force_close_closes_window(qtbot):
    view = TaviView()
    qtbot.addWidget(view)
    view.show()

    view.force_close()

    assert not view.isVisible()


def test_force_close_sets_force_closing_flag(qtbot):
    view = TaviView()
    qtbot.addWidget(view)
    view.show()

    view.force_close()

    assert view._force_closing is True


# ---------------------------------------------------------------------------
# TaviView — exit_message_box
# ---------------------------------------------------------------------------


def test_exit_message_box_returns_false_when_user_clicks_no(qtbot, monkeypatch):
    view = TaviView()
    qtbot.addWidget(view)

    monkeypatch.setattr(QMessageBox, "question", lambda *args, **kwargs: QMessageBox.No)
    assert view.exit_message_box() is False


def test_exit_message_box_returns_true_when_user_clicks_yes(qtbot, monkeypatch):
    view = TaviView()
    qtbot.addWidget(view)

    monkeypatch.setattr(QMessageBox, "question", lambda *args, **kwargs: QMessageBox.Yes)
    assert view.exit_message_box() is True
