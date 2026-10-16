"""Tests for FitEntry and FitMember."""

from tavi.library.data.fit_entry import FitEntry, FitMember, ParamField, PeakField
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
