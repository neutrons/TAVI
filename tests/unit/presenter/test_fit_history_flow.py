"""End-to-end tests for fit undo/redo across project tree selections, with real models and presenters."""

from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from tavi.backend.model.fit_model import FitModel
from tavi.backend.model.plot_model import PlotModel
from tavi.backend.model.tavi_project_model import TaviProjectModel
from tavi.frontend.presenter.fitting_presenter import FittingPresenter
from tavi.frontend.presenter.plotter_presenter import PlotterPresenter
from tavi.library.data.scan import UUID, Provenance, RawScan, ScanData, ScanMetadata, TaviMetadata
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.presenter_event import FocusEvent


def make_scan(uuid_val) -> RawScan:
    x = np.linspace(-5, 5, 101)
    y = 5.0 * np.exp(-(x**2) / 2)
    return RawScan(
        uuid=UUID(value=uuid_val),
        data=ScanData(data={"qh": x.tolist(), "detector": y.tolist()}),
        metadata=ScanMetadata(),
        tavimeta=TaviMetadata(default_axis=("qh", "detector"), friendly_name=uuid_val, friendly_path="/exp1"),
        prov=Provenance(raw_file="scan.dat", contributing_scans={UUID(value=uuid_val): 1}),
    )


@pytest.fixture
def app(qtbot):
    """Wire the project, plot and fit models to the plotter and fitting presenters, the way MainPresenter does."""
    with patch("tavi.backend.model.tavi_project_model.RawScanLoadController"):
        with patch("tavi.backend.model.tavi_project_model.Resource"):
            project = TaviProjectModel(MagicMock())
    for uuid_val in ("scan-001", "scan-002"):
        scan = make_scan(uuid_val)
        project.tavi_data.raw_scans[scan.uuid] = scan
    data = project.tavi_data
    plot_model = PlotModel(data.plots, data.raw_scans)
    fit_model = FitModel(data.raw_scans, data.fits)
    plotter = PlotterPresenter(plot_model)
    fitting = FittingPresenter(fit_model)
    qtbot.addWidget(plotter._view)
    qtbot.addWidget(fitting._view)
    return project, fitting._view


def select(*uuid_vals):
    EventBroker().publish(FocusEvent(ids=[UUID(value=v) for v in uuid_vals]))


def center_row(view):
    return view.peak_panels[0].peak_table.rows[1]


def fit_with_range_max(view, range_max):
    view.max_edit.setText(range_max)
    view.perform_fit_btn.click()


def history_buttons(view):
    return view.undo_fit_btn.isEnabled(), view.redo_fit_btn.isEnabled()


def test_refit_then_undo_restores_the_earlier_member(app):
    project, view = app
    select("scan-001")
    fit_with_range_max(view, "5")
    fit_with_range_max(view, "1")
    assert history_buttons(view) == (True, False)

    view.undo_fit_btn.click()

    assert view.max_edit.text() == "5"
    assert history_buttons(view) == (False, True)


def test_selecting_another_raw_scan_resets_fields_and_disables_history(app):
    project, view = app
    select("scan-001")
    fit_with_range_max(view, "5")
    center_row(view).fix_check.setChecked(True)
    fit_with_range_max(view, "1")

    select("scan-002")

    assert history_buttons(view) == (False, False)
    assert center_row(view).value_edit.text() == "1"
    assert not center_row(view).fix_check.isChecked()
    assert view.chi2_edit.text() == ""


def test_reselecting_a_raw_scan_starts_a_new_fit_with_its_own_history(app):
    project, view = app
    select("scan-001")
    fit_with_range_max(view, "5")
    fit_with_range_max(view, "1")

    select("scan-002")
    select("scan-001")
    assert history_buttons(view) == (False, False)
    fit_with_range_max(view, "3")

    assert len(project.tavi_data.fits) == 2
    assert history_buttons(view) == (False, False)


def test_selecting_a_fit_from_the_tree_brings_back_its_history(app):
    project, view = app
    select("scan-001")
    fit_with_range_max(view, "5")
    fit_with_range_max(view, "1")
    (fit_uuid,) = project.tavi_data.fits

    select("scan-002")
    select(fit_uuid.value)

    assert history_buttons(view) == (True, False)
    view.undo_fit_btn.click()
    assert view.max_edit.text() == "5"
