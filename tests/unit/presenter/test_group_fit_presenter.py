"""Tests for GroupFitPresenter."""

from unittest.mock import MagicMock

import numpy as np
import pytest

from tavi.backend.model.group_fit_model import GroupFitModel
from tavi.frontend.presenter.group_fit_presenter import GroupFitPresenter
from tavi.frontend.view.group_fit_view import GroupFitWindow
from tavi.library.data.fit_entry import FitCurve, FitData, FitSpec, ParamField, PeakField
from tavi.library.data.group_fit_entry import GroupFitOutcome, GroupFitRequest
from tavi.library.data.scan import UUID, Provenance, RawScan, ScanData, ScanMetadata, TaviMetadata
from tavi.meta.event.event_broker import EventBroker
from tavi.meta.event.type.model_event import RemoveProcessedScanEvent, RemoveRawScanEvent
from tavi.meta.event.type.presenter_event import (
    GroupScansEvent,
    SaveFitSummaryEvent,
    SyncGroupBackgroundParamsEvent,
    SyncGroupCosineParamsEvent,
    SyncGroupFitEvent,
    SyncGroupPeakParamsEvent,
)


def make_param(value="") -> ParamField:
    return ParamField(value=str(value), fixed=False, minimum="", maximum="")


def make_spec(center="") -> FitSpec:
    return FitSpec(
        range_min="-5",
        range_max="5",
        background="None",
        background_constant=make_param(0),
        peaks=[PeakField(shape="Gaussian", amplitude=make_param(), center=make_param(center), fwhm=make_param())],
    )


def make_scans(*names) -> list[tuple[UUID, str]]:
    return [(UUID(value=name), name) for name in names]


def make_raw_scan(name, center) -> RawScan:
    x = np.linspace(-5, 5, 201)
    detector = 50.0 * np.exp(-((x - center) ** 2) / (2 * 0.5**2)) + 1.0
    return RawScan(
        uuid=UUID(value=name),
        data=ScanData(data={"qh": x.tolist(), "detector": detector.tolist()}),
        metadata=ScanMetadata(),
        tavimeta=TaviMetadata(default_axis=("qh", "detector"), friendly_name=name, friendly_path="/exp1"),
        prov=Provenance(raw_file=f"{name}.dat", contributing_scans={UUID(value=name): 1}),
    )


def make_outcome(uuid) -> GroupFitOutcome:
    return GroupFitOutcome(
        source_scan_uuid=uuid,
        spec=make_spec(),
        data=FitData(x=[0.0, 1.0], y=[1.0, 2.0], err=[0.1, 0.1]),
        curve=FitCurve(source_scan_uuid=uuid, scan_name=uuid.value, x=[0.0, 1.0], best_fit=[1.0, 2.0]),
    )


def group(*names):
    EventBroker().publish(GroupScansEvent(scans=make_scans(*names)))


@pytest.fixture
def presenter(qtbot):
    p = GroupFitPresenter(MagicMock())
    qtbot.addWidget(p.view())
    return p


def last_request(presenter) -> GroupFitRequest:
    return presenter._model.perform_group_fit.call_args[0][0]


# ---------------------------------------------------------------------------
# Wiring
# ---------------------------------------------------------------------------


def test_view_is_a_group_fit_window(presenter):
    assert isinstance(presenter.view(), GroupFitWindow)


@pytest.mark.parametrize(
    ("event_type", "handler"),
    [
        (GroupScansEvent, "handle_group_scans"),
        (SyncGroupFitEvent, "handle_sync_group_fit"),
        (RemoveRawScanEvent, "handle_scan_removed"),
        (RemoveProcessedScanEvent, "handle_scan_removed"),
        (SyncGroupPeakParamsEvent, "handle_peak_params_suggested"),
        (SyncGroupCosineParamsEvent, "handle_cosine_params_suggested"),
        (SyncGroupBackgroundParamsEvent, "handle_background_params_suggested"),
    ],
)
def test_registers_for_its_events(presenter, event_type, handler):
    assert getattr(presenter, handler) in EventBroker().registry[event_type]


# ---------------------------------------------------------------------------
# Scans added from the project tree
# ---------------------------------------------------------------------------


def test_group_scans_event_lists_the_scans_and_opens_the_window(presenter):
    group("s1", "s2")

    assert presenter.view().isVisible()
    assert list(presenter.view().fit_specs()) == [UUID(value="s1"), UUID(value="s2")]


def test_group_scans_events_accumulate(presenter):
    group("s1")
    group("s2")

    assert list(presenter.view().fit_specs()) == [UUID(value="s1"), UUID(value="s2")]


@pytest.mark.parametrize("event_type", [RemoveRawScanEvent, RemoveProcessedScanEvent])
def test_a_scan_leaving_the_project_leaves_the_group(presenter, event_type):
    group("s1", "s2")

    EventBroker().publish(event_type(uuid=UUID(value="s1")))

    assert list(presenter.view().fit_specs()) == [UUID(value="s2")]


# ---------------------------------------------------------------------------
# Group Fit button
# ---------------------------------------------------------------------------


def test_group_fit_sends_every_scans_spec_and_the_columns(presenter):
    group("s1", "s2")
    view = presenter.view()
    view.set_columns("qh", "detector")
    view.fitting_view.min_edit.setText("-3")

    view.fitting_view.perform_fit_btn.click()

    request = last_request(presenter)
    assert list(request.specs) == [UUID(value="s1"), UUID(value="s2")]
    assert request.specs[UUID(value="s1")].range_min == "-3"
    assert (request.x_col, request.y_col) == ("qh", "detector")


def test_a_blank_x_column_is_sent_as_each_scans_default_axis(presenter):
    group("s1")
    presenter.view().set_columns("", "detector")

    presenter.view().fitting_view.perform_fit_btn.click()

    assert last_request(presenter).x_col is None


# ---------------------------------------------------------------------------
# Outcomes
# ---------------------------------------------------------------------------


def test_sync_group_fit_draws_the_selected_scans(presenter):
    group("s1", "s2")

    EventBroker().publish(
        SyncGroupFitEvent(x_col=None, y_col="detector", outcomes=[make_outcome(UUID(value=n)) for n in ("s1", "s2")])
    )

    assert presenter.view().axes.get_legend_handles_labels()[1] == ["s1 fit", "s1"]


def test_group_fit_end_to_end_overplots_the_fitted_scans(qtbot):
    """Tree -> window -> Group Fit -> GroupFitModel -> SyncGroupFitEvent -> overplot, with the real model."""
    raw_scans = {UUID(value=name): make_raw_scan(name, center) for name, center in (("s1", -1.0), ("s2", 0.0))}
    presenter = GroupFitPresenter(GroupFitModel(raw_scans, {}))
    qtbot.addWidget(presenter.view())
    view = presenter.view()
    group("s1", "s2")
    view.set_columns("", "detector")
    view.fitting_view.peak_panels[0].peak_table.rows[1].value_edit.setText("-1")

    view.fitting_view.perform_fit_btn.click()
    view.scan_list.selectAll()

    assert sorted(view.axes.get_legend_handles_labels()[1]) == ["s1", "s1 fit", "s2", "s2 fit"]
    assert view.axes.get_xlabel() == "default x axis"


# ---------------------------------------------------------------------------
# Summarize fit results
# ---------------------------------------------------------------------------


def test_generate_before_a_group_fit_asks_nothing(presenter):
    group("s1")

    presenter.view().generate_summary_btn.click()

    presenter._model.summarize_group_fit.assert_not_called()


def test_generate_sends_the_windows_summary_request(presenter):
    group("s1")
    view = presenter.view()
    view.summary_param_combo.addItem("peak1 center")
    view.summary_x_edit.setText("temp")
    view.generate_summary_btn.setEnabled(True)

    view.generate_summary_btn.click()

    request = presenter._model.summarize_group_fit.call_args[0][0]
    assert (request.parameter, request.x_col) == ("peak1 center", "temp")


def test_generate_end_to_end_saves_the_fitted_centers(qtbot):
    """Group Fit -> pick peak1 center vs qh -> Generate -> SaveFitSummaryEvent, with the real model."""
    raw_scans = {UUID(value=name): make_raw_scan(name, center) for name, center in (("s1", -1.0), ("s2", 0.0))}
    presenter = GroupFitPresenter(GroupFitModel(raw_scans, {}))
    qtbot.addWidget(presenter.view())
    view = presenter.view()
    saved = []
    EventBroker().register(SaveFitSummaryEvent, saved.append)
    group("s1", "s2")
    view.set_columns("", "detector")
    view.fitting_view.peak_panels[0].peak_table.rows[1].value_edit.setText("-1")
    view.fitting_view.perform_fit_btn.click()

    view.summary_param_combo.setCurrentText("peak1 center")
    view.summary_x_edit.setText("qh")
    view.generate_summary_btn.click()

    assert len(saved) == 1
    points = saved[0].summary.points
    assert [point.scan_name for point in points] == ["s1", "s2"]
    assert [point.y for point in points] == pytest.approx([-1.0, 0.0], abs=1e-3)
    assert [point.x for point in points] == pytest.approx([0.0, 0.0], abs=1e-9)


# ---------------------------------------------------------------------------
# Suggest Params.
# ---------------------------------------------------------------------------


def test_peak_suggest_asks_the_model_about_the_scan_on_screen(presenter):
    group("s1", "s2")
    view = presenter.view()
    view.set_columns("", "detector")
    view.scan_list.setCurrentRow(1)

    view.fitting_view.peak_panels[0].suggest_btn.click()

    request = presenter._model.suggest_peak_params.call_args[0][0]
    assert (request.source_scan_uuid, request.y_col, request.shape) == (UUID(value="s2"), "detector", "Gaussian")


def test_background_suggest_asks_the_model(presenter):
    group("s1")

    presenter.view().fitting_view.suggest_background_btn.click()

    assert presenter._model.suggest_background_params.call_args[0][0].background == "Linear"


def test_suggest_without_a_scan_asks_nothing(presenter):
    presenter.handle_suggest_params_clicked()
    presenter.handle_suggest_background_clicked()

    presenter._model.suggest_peak_params.assert_not_called()
    presenter._model.suggest_background_params.assert_not_called()


def test_guess_events_fill_the_window(presenter):
    group("s1")
    uuid = UUID(value="s1")

    EventBroker().publish(SyncGroupPeakParamsEvent(source_scan_uuid=uuid, amplitude=9.0, center=0.5, fwhm=0.2))
    EventBroker().publish(SyncGroupBackgroundParamsEvent(source_scan_uuid=uuid, slope=0.1, intercept=2.0))

    view = presenter.view().fitting_view
    assert [row.value_edit.text() for row in view.peak_panels[0].peak_table.rows] == ["9", "0.5", "0.2"]
    assert [row.value_edit.text() for row in view.background_table.rows] == ["0.1", "2"]


def test_suggest_end_to_end_fills_the_panel_from_that_scans_data(qtbot):
    raw_scans = {UUID(value=name): make_raw_scan(name, center) for name, center in (("s1", -1.0), ("s2", 1.5))}
    presenter = GroupFitPresenter(GroupFitModel(raw_scans, {}))
    qtbot.addWidget(presenter.view())
    view = presenter.view()
    group("s1", "s2")
    view.set_columns("", "detector")
    view.scan_list.setCurrentRow(1)

    view.fitting_view.peak_panels[0].suggest_btn.click()

    center = float(view.fitting_view.peak_panels[0].peak_table.rows[1].value_edit.text())
    assert center == pytest.approx(1.5, abs=0.1)
