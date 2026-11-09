"""Tests for CombineOps, the one entry point that appends or rebins and stores the result."""

import pytest

from tavi.library.data.algorithms.append_ops import AppendOp
from tavi.library.data.algorithms.bin_by_tol import ToleranceRebinOps
from tavi.library.data.algorithms.equal_step_bin import EqualStepRebinOps
from tavi.library.data.combine import CombineOps
from tavi.library.data.enum.rebin_mode import RebinMode
from tavi.library.data.scan import UUID, ProcessedScan, Provenance, RawScan, ScanData, ScanMetadata, TaviMetadata
from tavi.library.data.tavi_data import TaviData

COLUMNS = ["qh", "detector", "monitor"]


def make_raw_scan(
    uuid_val="scan-001",
    name="HB1A_exp0004_scan0001",
    data=None,
    default_axis=("qh", "detector"),
) -> RawScan:
    """Build a RawScan standing in for one loaded from disk."""
    return RawScan(
        uuid=UUID(value=uuid_val),
        data=ScanData(data=data if data is not None else {"qh": [1.0], "detector": [100.0], "monitor": [1000.0]}),
        metadata=ScanMetadata(),
        tavimeta=TaviMetadata(
            default_axis=default_axis, friendly_name=name, friendly_path="/IPTS-1", normalization=("monitor", 1.0)
        ),
        prov=Provenance(raw_file=f"{name}.dat", contributing_scans={UUID(value=uuid_val): 1}),
    )


def make_countable_scan(uuid_val="scan-001", data=None) -> RawScan:
    """A scan carrying the columns the rebin paths ask for: an x, raw counts and a monitor."""
    return make_raw_scan(
        uuid_val,
        data=data
        if data is not None
        else {"qh": [1.0, 1.4, 2.1], "detector": [10.0, 20.0, 30.0], "monitor": [1000.0] * 3},
    )


def make_tavi_data(*scans) -> TaviData:
    """Put scans in the pool CombineOps resolves uuids against."""
    return TaviData(raw_scans={scan.uuid: scan for scan in scans})


def make_processed_scan(tavi_data, scan) -> ProcessedScan:
    """Append a scan to itself and register the result, standing in for any already-reduced origin."""
    processed = AppendOp(tavi_data, [scan.uuid], COLUMNS).exec()
    tavi_data.processed_scans[processed.uuid] = processed
    return processed


# --- which operation the parameters resolve to ------------------------------------------


def test_combine_with_a_tolerance_runs_the_tolerance_rebin():
    """A tol names the tolerance rebin, and combine must hand back exactly what it builds."""
    scan = make_countable_scan()
    direct = ToleranceRebinOps(make_tavi_data(scan), [scan.uuid], COLUMNS, norm="monitor", tol=0.5).exec()

    combined = CombineOps(make_tavi_data(scan), [scan.uuid], COLUMNS, norm="monitor", tol=0.5).exec()

    assert combined.data.data == direct.data.data
    assert combined.tavimeta.friendly_name == direct.tavimeta.friendly_name


def test_combine_with_a_step_runs_the_equal_step_rebin():
    """A step names the equal-step rebin, whose x values sit on the grid rather than on the data."""
    scan = make_countable_scan()
    direct = EqualStepRebinOps(make_tavi_data(scan), [scan.uuid], COLUMNS, norm="monitor", step=1.0).exec()

    combined = CombineOps(make_tavi_data(scan), [scan.uuid], COLUMNS, norm="monitor", step=1.0).exec()

    assert combined.data.data == direct.data.data
    assert combined.tavimeta.friendly_name == direct.tavimeta.friendly_name


def test_combine_without_a_tolerance_or_a_step_appends():
    """Neither parameter given is a plain append - every row survives and no error column appears."""
    first = make_countable_scan("scan-001", data={"qh": [1.0], "detector": [10.0], "monitor": [1000.0]})
    second = make_countable_scan("scan-002", data={"qh": [1.05], "detector": [20.0], "monitor": [1000.0]})

    combined = CombineOps(make_tavi_data(first, second), [first.uuid, second.uuid], COLUMNS).exec()

    assert combined.data.qh == [1.0, 1.05]
    assert list(combined.data.data) == COLUMNS


def test_combine_passes_the_grid_and_the_scale_through_to_the_equal_step_rebin():
    """start, finish and multiplication are the equal-step rebin's, and must arrive unaltered."""
    scan = make_countable_scan()
    direct = EqualStepRebinOps(
        make_tavi_data(scan), [scan.uuid], COLUMNS, norm="monitor", step=1.0, start=1.0, finish=3.0, multiplication=1e4
    ).exec()

    combined = CombineOps(
        make_tavi_data(scan), [scan.uuid], COLUMNS, norm="monitor", step=1.0, start=1.0, finish=3.0, multiplication=1e4
    ).exec()

    assert combined.data.data == direct.data.data
    assert combined.tavimeta.normalization == ("monitor", 1e4)


def test_combine_merges_two_scans_across_their_seam():
    """Binning the combination rather than each origin is what a combine of repeated passes is for."""
    first = make_raw_scan("scan-001", data={"qh": [1.0], "detector": [100.0], "monitor": [1000.0]})
    second = make_raw_scan("scan-002", data={"qh": [1.2], "detector": [140.0], "monitor": [1000.0]})
    tavi_data = make_tavi_data(first, second)

    combined = CombineOps(
        tavi_data, [first.uuid, second.uuid], COLUMNS, norm="monitor", step=1.0, multiplication=1000
    ).exec()

    assert combined.data.detector == pytest.approx([1000 * 240 / 2000])
    assert combined.prov.contributing_scans == {first.uuid: 1, second.uuid: 1}


def test_combine_accepts_a_single_scan():
    """One origin is a combination of one, which the rebins are as happy with as the append."""
    scan = make_countable_scan()

    combined = CombineOps(make_tavi_data(scan), [scan.uuid], COLUMNS, norm="monitor", tol=0.5).exec()

    assert combined.prov.contributing_scans == {scan.uuid: 1}


# --- storing ----------------------------------------------------------------------------


def test_exec_registers_the_result_in_processed_scans():
    """Combining is how a derived scan comes to exist in a project, so exec stores as well as returns."""
    scan = make_countable_scan()
    tavi_data = make_tavi_data(scan)

    combined = CombineOps(tavi_data, [scan.uuid], COLUMNS, norm="monitor", tol=0.5).exec()

    assert tavi_data.processed_scans[combined.uuid] is combined
    assert tavi_data.fetch_by_uuid(combined.uuid) is combined


def test_exec_leaves_the_raw_scans_alone():
    """Storing the result must not touch the pool it was derived from."""
    scan = make_countable_scan()
    tavi_data = make_tavi_data(scan)

    CombineOps(tavi_data, [scan.uuid], COLUMNS, norm="monitor", tol=0.5).exec()

    assert set(tavi_data.raw_scans) == {scan.uuid}
    assert scan.data.qh == [1.0, 1.4, 2.1]


def test_two_combines_of_the_same_scans_are_stored_side_by_side():
    """Each result carries a fresh uuid, so a second normalization of one scan never overwrites the first."""
    scan = make_countable_scan()
    tavi_data = make_tavi_data(scan)

    once = CombineOps(tavi_data, [scan.uuid], COLUMNS, norm="monitor", tol=0.5).exec()
    twice = CombineOps(tavi_data, [scan.uuid], COLUMNS, norm="monitor", step=1.0).exec()

    assert once.uuid != twice.uuid
    assert set(tavi_data.processed_scans) == {once.uuid, twice.uuid}


def test_validate_on_its_own_stores_nothing():
    """A GUI may offer the operation before the user commits to it, which must leave no trace."""
    scan = make_countable_scan()
    tavi_data = make_tavi_data(scan)

    CombineOps(tavi_data, [scan.uuid], COLUMNS, norm="monitor", tol=0.5).validate()

    assert tavi_data.processed_scans == {}


# --- inspection -------------------------------------------------------------------------


@pytest.mark.parametrize(
    "kwargs, expected",
    [
        ({"norm": "monitor", "tol": 0.5}, RebinMode.TOLERANCE),
        ({"norm": "monitor", "step": 1.0}, RebinMode.EQUAL_STEP),
        ({}, RebinMode.NONE),
    ],
)
def test_mode_reports_what_the_parameters_resolve_to(kwargs, expected):
    """A caller can ask a built op what it will do before running it."""
    scan = make_countable_scan()

    assert CombineOps(make_tavi_data(scan), [scan.uuid], COLUMNS, **kwargs).mode is expected


def test_mode_and_exec_still_agree_after_a_parameter_is_reassigned():
    """The delegate is built afresh at exec, so mode never reports one combination while another runs."""
    scan = make_countable_scan()
    tavi_data = make_tavi_data(scan)
    op = CombineOps(tavi_data, [scan.uuid], COLUMNS, norm="monitor", step=1.0)
    op.validate()

    op.step, op.tol = None, 0.5

    assert op.mode is RebinMode.TOLERANCE
    assert "_tol0.5" in op.exec().tavimeta.friendly_name


def test_exec_picks_up_a_reassigned_value_as_the_delegates_do():
    """A combine reads its parameters when it runs, matching the ops it delegates to."""
    scan = make_countable_scan()
    op = CombineOps(make_tavi_data(scan), [scan.uuid], COLUMNS, norm="monitor", step=1.0)
    op.validate()

    op.step = 0.5

    assert "_rebin0.5" in op.exec().tavimeta.friendly_name


def test_columns_and_size_report_what_the_op_was_built_with():
    scan = make_countable_scan()

    op = CombineOps(make_tavi_data(scan), [scan.uuid], COLUMNS, norm="monitor", tol=0.5)

    assert list(op.columns) == COLUMNS
    assert op.size == 1


# --- what combine refuses ---------------------------------------------------------------


def test_combine_refuses_a_tolerance_and_a_step_together():
    """Two mutually exclusive grids, with nothing to say which one was meant."""
    scan = make_countable_scan()

    with pytest.raises(ValueError, match="either a tolerance or a step"):
        CombineOps(make_tavi_data(scan), [scan.uuid], COLUMNS, norm="monitor", tol=0.5, step=1.0).validate()


def test_combine_refuses_no_scans_at_all():
    """An append of nothing yields empty columns rather than failing, which is not a combination."""
    with pytest.raises(ValueError, match="at least one scan"):
        CombineOps(TaviData(), [], COLUMNS).validate()


@pytest.mark.parametrize(
    "kwargs",
    [{"norm": "monitor", "tol": 0.5}, {"norm": "monitor", "step": 1.0}, {}],
)
def test_combine_refuses_an_already_processed_origin_in_every_mode(kwargs):
    """Raw counts against a monitor are meaningless on data something has already reduced."""
    scan = make_countable_scan()
    tavi_data = make_tavi_data(scan)
    processed = make_processed_scan(tavi_data, scan)

    with pytest.raises(ValueError, match="only operate on raw data"):
        CombineOps(tavi_data, [processed.uuid], COLUMNS, **kwargs).validate()


@pytest.mark.parametrize("edge", ["start", "finish"])
def test_combine_refuses_a_grid_edge_without_a_step(edge):
    """An edge is a position on a grid, and only a step lays one down."""
    scan = make_countable_scan()

    with pytest.raises(ValueError, match="only applies to an equal-step rebin"):
        CombineOps(make_tavi_data(scan), [scan.uuid], COLUMNS, norm="monitor", tol=0.5, **{edge: 1.0}).validate()


def test_combine_refuses_a_normalization_it_has_no_rebin_to_apply():
    """Dropping it silently would hand back an unnormalized scan to a caller who asked for a normalized one."""
    scan = make_countable_scan()

    with pytest.raises(ValueError, match="needs a rebin to weigh by it"):
        CombineOps(make_tavi_data(scan), [scan.uuid], COLUMNS, norm="monitor").validate()


def test_combine_refuses_a_multiplication_it_has_no_rebin_to_apply():
    scan = make_countable_scan()

    with pytest.raises(ValueError, match="needs a rebin to scale"):
        CombineOps(make_tavi_data(scan), [scan.uuid], COLUMNS, multiplication=1000).validate()


def test_combine_still_defers_to_the_delegate_for_its_own_checks():
    """The parameters belonging to the mode is combine's business; whether they are usable is the op's."""
    scan = make_countable_scan()

    with pytest.raises(ValueError, match="must be one of the requested columns"):
        CombineOps(make_tavi_data(scan), [scan.uuid], COLUMNS, norm="not_a_column", tol=0.5).validate()


def test_combine_refuses_a_scan_that_is_not_loaded():
    scan = make_countable_scan()

    with pytest.raises(KeyError):
        CombineOps(make_tavi_data(scan), [UUID(value="not-loaded")], COLUMNS, norm="monitor", tol=0.5).validate()


def test_a_failed_combine_stores_nothing():
    """exec validates first, so a refused combination must leave the pool as it found it."""
    scan = make_countable_scan()
    tavi_data = make_tavi_data(scan)

    with pytest.raises(ValueError):
        CombineOps(tavi_data, [scan.uuid], COLUMNS, norm="monitor", tol=0.5, step=1.0).exec()

    assert tavi_data.processed_scans == {}
