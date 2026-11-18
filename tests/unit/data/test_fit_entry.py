"""Tests for FitEntry, FitMember and FitResultSummary."""

from tavi.library.data.fit_entry import (
    FitEntry,
    FitMember,
    FitResultSummary,
    FitSummary,
    FitSummaryPoint,
    ParamField,
    PeakField,
    PeakResult,
)
from tavi.library.data.plot import PlotSeries
from tavi.library.data.scan import UUID


def make_param(value=0) -> ParamField:
    return ParamField(value=str(value), fixed=False, minimum="", maximum="")


def make_series(uuid_val="scan-001") -> PlotSeries:
    return PlotSeries(
        source_scan_uuid=UUID(value=uuid_val),
        scan_name=uuid_val,
        normalized_by=None,
        x_name="qh",
        y_name="en",
        error_name="error",
    )


def make_member(uuid_val="scan-001", range_max="10") -> FitMember:
    return FitMember(
        series=make_series(uuid_val),
        range_min="0",
        range_max=range_max,
        background="None",
        background_constant=make_param(0),
        peaks=[PeakField(shape="Gaussian", amplitude=make_param(1), center=make_param(0), fwhm=make_param(1))],
    )


def make_fit_entry(*scan_uuids: str) -> FitEntry:
    return FitEntry(uuid=UUID(value="fit-001"), members=[make_member(u) for u in scan_uuids])


def sources(fit: FitEntry) -> list[str]:
    return [member.source_scan_uuid.value for member in fit.members]


def test_member_source_scan_uuid_is_its_series_source():
    assert make_member("scan-007").source_scan_uuid == UUID(value="scan-007")


def test_member_for_finds_the_member_for_a_scan():
    fit = make_fit_entry("scan-001", "scan-002")

    assert fit.member_for(UUID(value="scan-002")) is fit.members[1]


def test_member_for_returns_none_for_an_uncovered_scan():
    assert make_fit_entry("scan-001").member_for(UUID(value="scan-other")) is None


def test_with_members_replaces_the_member_for_the_same_scan_in_place():
    fit = make_fit_entry("scan-001", "scan-002", "scan-003")

    merged = fit.with_members([make_member("scan-002", range_max="5")])

    assert sources(merged) == ["scan-001", "scan-002", "scan-003"]
    assert [m.range_max for m in merged.members] == ["10", "5", "10"]


def test_with_members_appends_a_member_for_a_new_scan():
    fit = make_fit_entry("scan-001")

    merged = fit.with_members([make_member("scan-002")])

    assert sources(merged) == ["scan-001", "scan-002"]


def test_with_members_returns_a_copy_and_keeps_identity():
    fit = make_fit_entry("scan-001")
    fit = fit.model_copy(update={"name": "run_Fit"})

    merged = fit.with_members([make_member("scan-001", range_max="5")])

    assert merged is not fit
    assert fit.members[0].range_max == "10"
    assert (merged.uuid, merged.name) == (fit.uuid, fit.name)


def test_without_scans_returns_self_when_nothing_matches():
    fit = make_fit_entry("scan-001", "scan-002")

    assert fit.without_scans({UUID(value="scan-other")}) is fit


def test_without_scans_drops_only_matching_members():
    fit = make_fit_entry("scan-001", "scan-002", "scan-003")

    pruned = fit.without_scans({UUID(value="scan-002")})

    assert sources(pruned) == ["scan-001", "scan-003"]
    assert pruned.uuid == fit.uuid
    assert sources(fit) == ["scan-001", "scan-002", "scan-003"]


def test_without_scans_returns_none_once_no_member_is_left():
    fit = make_fit_entry("scan-001", "scan-002")

    assert fit.without_scans({UUID(value="scan-001"), UUID(value="scan-002")}) is None


def make_peak_result(amplitude=10.0, center=1.0, fwhm=0.5) -> PeakResult:
    return PeakResult(amplitude=amplitude, amplitude_err=0.1, center=center, center_err=0.01, fwhm=fwhm, fwhm_err=0.02)


def make_cosine_result() -> PeakResult:
    return PeakResult(
        amplitude=2.0,
        amplitude_err=0.2,
        frequency=3.0,
        frequency_err=0.3,
        phase=0.4,
        phase_err=0.04,
        offset=5.0,
        offset_err=0.5,
    )


def test_parameters_names_a_peaks_amplitude_its_integrated_intensity():
    summary = FitResultSummary(reduced_chi_squared=1.2, peaks=[make_peak_result()])

    assert summary.parameters() == {
        "peak1 Int. Inten.": (10.0, 0.1),
        "peak1 center": (1.0, 0.01),
        "peak1 FWHM": (0.5, 0.02),
        "reduced χ²": (1.2, None),
    }


def test_parameters_number_every_peak_and_add_the_background():
    summary = FitResultSummary(
        reduced_chi_squared=1.0,
        peaks=[make_peak_result(), make_peak_result(center=2.0)],
        background_constant=3.0,
        background_constant_err=0.3,
        background_slope=0.5,
        background_slope_err=0.05,
    )

    assert list(summary.parameters()) == [
        "peak1 Int. Inten.",
        "peak1 center",
        "peak1 FWHM",
        "peak2 Int. Inten.",
        "peak2 center",
        "peak2 FWHM",
        "bg slope",
        "bg constant",
        "reduced χ²",
    ]
    assert summary.parameters()["peak2 center"] == (2.0, 0.01)
    assert summary.parameters()["bg constant"] == (3.0, 0.3)


def test_parameters_keep_a_cosines_own_names():
    summary = FitResultSummary(reduced_chi_squared=1.0, peaks=[make_cosine_result()])

    assert list(summary.parameters()) == [
        "peak1 amplitude",
        "peak1 frequency",
        "peak1 phase",
        "peak1 offset",
        "reduced χ²",
    ]


def make_fit_summary(uuid_val="summary-001") -> FitSummary:
    points = [
        FitSummaryPoint(source_scan_uuid=UUID(value="s0"), scan_name="s0", x=10.0, x_err=0.5, y=1.0, y_err=0.1),
        FitSummaryPoint(source_scan_uuid=UUID(value="s1"), scan_name="s1", x=20.0, x_err=0.4, y=1.5, y_err=None),
    ]
    return FitSummary(
        uuid=UUID(value=uuid_val), name="peak1 center vs temp", x_label="temp", y_label="peak1 center", points=points
    )


def test_fit_summary_columns_list_each_scans_point_with_its_uncertainties():
    assert make_fit_summary().columns() == {
        "scan": ["s0", "s1"],
        "temp": [10.0, 20.0],
        "temp err": [0.5, 0.4],
        "peak1 center": [1.0, 1.5],
        "peak1 center err": [0.1, None],
    }
