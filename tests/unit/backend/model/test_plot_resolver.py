import numpy.testing as npt
import pytest

from tavi.backend.model.plot_resolver import (
    channel_pairs,
    expand_series,
    find_series_by_source,
    matching_columns,
    resolve_series,
    scans_for_plots,
)
from tavi.library.data.plot import Plot, PlotSeries
from tavi.library.data.scan import UUID, Provenance, RawScan, ScanData, ScanMetadata, TaviMetadata


def make_scan(uuid_val="scan-001", data=None) -> RawScan:
    if data is None:
        data = {"qh": [1.0, 2.0, 3.0], "en": [4.0, 9.0, 16.0]}
    return RawScan(
        uuid=UUID(value=uuid_val),
        data=ScanData(data=data),
        metadata=ScanMetadata(),
        tavimeta=TaviMetadata(default_axis=("qh", "en"), friendly_name="test_scan", friendly_path="/exp1"),
        prov=Provenance(raw_file="scan.dat", contributing_scans={UUID(value=uuid_val): 1}),
    )


def make_series(uuid_val="scan-001", x_name="qh", y_name="en", normalized_by=None, normalized_by_value=None) -> PlotSeries:
    return PlotSeries(
        source_scan_uuid=UUID(value=uuid_val),
        scan_name="test_scan",
        normalized_by=normalized_by,
        normalized_by_value=normalized_by_value,
        x_name=x_name,
        y_name=y_name,
        error_name="error",
    )


def test_resolve_series_returns_raw_x_and_y():
    scan = make_scan()
    series = make_series()

    x, y, _ = resolve_series(series, {scan.uuid: scan})

    npt.assert_array_equal(x, [1.0, 2.0, 3.0])
    npt.assert_array_equal(y, [4.0, 9.0, 16.0])


def test_resolve_series_error_is_sqrt_abs_y():
    scan = make_scan()
    series = make_series()

    _, y, err = resolve_series(series, {scan.uuid: scan})

    npt.assert_allclose(err, [2.0, 3.0, 4.0])


def test_resolve_series_unnormalized_when_normalized_by_is_none():
    scan = make_scan()
    series = make_series(normalized_by=None)

    _, y, _ = resolve_series(series, {scan.uuid: scan})

    npt.assert_array_equal(y, [4.0, 9.0, 16.0])


def test_resolve_series_divides_y_by_normalization_channel():
    scan = make_scan(data={"qh": [1.0, 2.0, 3.0], "en": [4.0, 9.0, 16.0], "monitor": [2.0, 3.0, 4.0]})
    series = make_series(normalized_by="monitor", normalized_by_value=1.0)

    _, y, _ = resolve_series(series, {scan.uuid: scan})

    npt.assert_allclose(y, [4.0 / 2.0, 9.0 / 3.0, 16.0 / 4.0])


def test_resolve_series_scales_by_normalized_by_value():
    scan = make_scan(data={"qh": [1.0, 2.0, 3.0], "en": [4.0, 9.0, 16.0], "monitor": [2.0, 3.0, 4.0]})
    series = make_series(normalized_by="monitor", normalized_by_value=10.0)

    _, y, _ = resolve_series(series, {scan.uuid: scan})

    npt.assert_allclose(y, [4.0 * 10.0 / 2.0, 9.0 * 10.0 / 3.0, 16.0 * 10.0 / 4.0])


def test_resolve_series_defaults_normalized_by_value_to_one_when_none():
    scan = make_scan(data={"qh": [1.0, 2.0, 3.0], "en": [4.0, 9.0, 16.0], "monitor": [2.0, 3.0, 4.0]})
    series = make_series(normalized_by="monitor", normalized_by_value=None)

    _, y, _ = resolve_series(series, {scan.uuid: scan})

    npt.assert_allclose(y, [4.0 / 2.0, 9.0 / 3.0, 16.0 / 4.0])


def test_resolve_series_propagates_err_through_normalization():
    scan = make_scan(data={"qh": [1.0, 2.0, 3.0], "en": [4.0, 9.0, 16.0], "monitor": [2.0, 6.0, 8.0]})
    series = make_series(normalized_by="monitor", normalized_by_value=10.0)

    _, _, err = resolve_series(series, {scan.uuid: scan})

    # The counting error is scaled the same way y is: sqrt(|y|) * weight / monitor.
    npt.assert_allclose(err, [10.0 * 2.0 / 2.0, 10.0 * 3.0 / 6.0, 10.0 * 4.0 / 8.0])


def test_scans_for_plots_returns_only_referenced_scans():
    scan_a = make_scan(uuid_val="scan-a")
    scan_b = make_scan(uuid_val="scan-b")
    plot = Plot(series=[make_series(uuid_val="scan-a")])

    result = scans_for_plots([plot], {scan_a.uuid: scan_a, scan_b.uuid: scan_b})

    assert set(result) == {scan_a.uuid}


def test_find_series_by_source_finds_a_series_within_a_multi_series_plot():
    """The point of this lookup: pick one series out of an otherwise-fused, multi-series plot."""
    plot = Plot(series=[make_series(uuid_val="scan-a"), make_series(uuid_val="scan-b")])

    result = find_series_by_source([plot], UUID(value="scan-b"))

    assert result is not None
    found_plot, found_series = result
    assert found_plot is plot
    assert found_series.source_scan_uuid == UUID(value="scan-b")


def test_find_series_by_source_searches_across_every_plot_in_the_batch():
    plot_a = Plot(series=[make_series(uuid_val="scan-a")])
    plot_b = Plot(series=[make_series(uuid_val="scan-b")])

    result = find_series_by_source([plot_a, plot_b], UUID(value="scan-b"))

    assert result is not None
    found_plot, _ = result
    assert found_plot is plot_b


def test_find_series_by_source_returns_none_when_no_series_matches():
    plot = Plot(series=[make_series(uuid_val="scan-a")])

    result = find_series_by_source([plot], UUID(value="not-focused"))

    assert result is None


def test_find_series_by_source_empty_plots_returns_none():
    assert find_series_by_source([], UUID(value="scan-a")) is None


# ---------------------------------------------------------------------------
# Wildcard channel names
# ---------------------------------------------------------------------------


def make_multichannel_scan(uuid_val="scan-001") -> RawScan:
    return make_scan(
        uuid_val=uuid_val,
        data={
            "qh": [1.0, 2.0],
            "detector_1": [2.0, 4.0],
            "detector_2": [6.0, 8.0],
            "detector_10": [10.0, 12.0],
            "monitor_1": [2.0, 2.0],
            "monitor_2": [2.0, 2.0],
            "monitor_10": [2.0, 2.0],
        },
    )


def test_matching_columns_without_a_wildcard_is_the_name_itself():
    assert matching_columns("detector_1", ["detector_1", "detector_2"]) == ["detector_1"]


def test_matching_columns_without_a_wildcard_matches_nothing_when_absent():
    assert matching_columns("detector_3", ["detector_1", "detector_2"]) == []


def test_matching_columns_takes_every_column_the_prefix_starts():
    assert matching_columns("detector*", ["detector_1", "monitor_1", "detector_2"]) == ["detector_1", "detector_2"]


def test_matching_columns_orders_numbered_channels_numerically():
    """detector_10 after detector_2, so the n-th y channel and n-th preset channel are the same detector."""
    assert matching_columns("detector*", ["detector_10", "detector_1", "detector_2"]) == [
        "detector_1",
        "detector_2",
        "detector_10",
    ]


def test_matching_columns_bare_wildcard_takes_everything():
    assert matching_columns("*", ["qh", "en"]) == ["en", "qh"]


def test_channel_pairs_without_normalization_pairs_each_match_with_none():
    pairs = channel_pairs(["detector_1", "detector_2", "qh"], "detector*", None)

    assert pairs == [("detector_1", None), ("detector_2", None)]


def test_channel_pairs_zips_two_wildcards_in_sorted_order():
    columns = ["detector_1", "detector_2", "monitor_1", "monitor_2"]

    pairs = channel_pairs(columns, "detector*", "monitor*")

    assert pairs == [("detector_1", "monitor_1"), ("detector_2", "monitor_2")]


def test_channel_pairs_reuses_a_plain_preset_channel_for_every_y_channel():
    columns = ["detector_1", "detector_2", "monitor"]

    pairs = channel_pairs(columns, "detector*", "monitor")

    assert pairs == [("detector_1", "monitor"), ("detector_2", "monitor")]


def test_channel_pairs_reuses_a_plain_y_channel_for_every_preset_channel():
    columns = ["detector", "monitor_1", "monitor_2"]

    pairs = channel_pairs(columns, "detector", "monitor*")

    assert pairs == [("detector", "monitor_1"), ("detector", "monitor_2")]


def test_channel_pairs_rejects_a_y_name_that_matches_nothing():
    with pytest.raises(ValueError, match="detector"):
        channel_pairs(["qh", "en"], "detector*", None)


def test_channel_pairs_rejects_a_preset_channel_that_matches_nothing():
    with pytest.raises(ValueError, match="monitor"):
        channel_pairs(["detector_1"], "detector*", "monitor*")


def test_channel_pairs_rejects_two_wildcards_of_different_lengths():
    columns = ["detector_1", "detector_2", "detector_3", "monitor_1", "monitor_2"]

    with pytest.raises(ValueError, match="paired up"):
        channel_pairs(columns, "detector*", "monitor*")


def test_channel_pairs_broadcasts_a_wildcard_that_matches_one_column():
    """One monitor for the whole bank is the normal case, and "monitor*" is a fine way to name it."""
    columns = ["detector_1", "detector_2", "monitor_1"]

    pairs = channel_pairs(columns, "detector*", "monitor*")

    assert pairs == [("detector_1", "monitor_1"), ("detector_2", "monitor_1")]


def test_expand_series_without_a_wildcard_returns_the_series_itself():
    scan = make_scan()
    series = make_series()

    assert expand_series(series, {scan.uuid: scan}) == [series]


def test_expand_series_makes_one_series_per_matched_channel():
    scan = make_multichannel_scan()
    series = make_series(y_name="detector*")

    expanded = expand_series(series, {scan.uuid: scan})

    assert [s.y_channel for s in expanded] == ["detector_1", "detector_2", "detector_10"]


def test_expand_series_pairs_each_y_channel_with_its_preset_channel():
    scan = make_multichannel_scan()
    series = make_series(y_name="detector*", normalized_by="monitor*")

    expanded = expand_series(series, {scan.uuid: scan})

    assert [(s.y_channel, s.normalized_by_channel) for s in expanded] == [
        ("detector_1", "monitor_1"),
        ("detector_2", "monitor_2"),
        ("detector_10", "monitor_10"),
    ]


def test_expand_series_keeps_the_names_as_typed_on_every_copy():
    """The wildcard is what the user asked for - it stays, so the fields and axis labels still show it."""
    scan = make_multichannel_scan()
    series = make_series(y_name="detector*", normalized_by="monitor*")

    expanded = expand_series(series, {scan.uuid: scan})

    assert {(s.y_name, s.normalized_by) for s in expanded} == {("detector*", "monitor*")}


def test_expand_series_leaves_an_unmatchable_name_alone():
    """Validation is update_fields' job; here an unmatched name must fail as a missing column does."""
    scan = make_scan()
    series = make_series(y_name="detector*")

    assert expand_series(series, {scan.uuid: scan}) == [series]


def test_expand_series_of_a_wildcard_matching_one_column_still_names_its_channel():
    scan = make_scan(data={"qh": [1.0], "detector_1": [2.0]})
    series = make_series(y_name="detector*")

    expanded = expand_series(series, {scan.uuid: scan})

    assert [s.channel_label for s in expanded] == ["detector_1"]


def test_resolve_series_of_an_expanded_copy_reads_its_own_channel():
    scan = make_multichannel_scan()
    series = make_series(y_name="detector*", normalized_by="monitor*", normalized_by_value=1.0)

    resolved = [resolve_series(s, {scan.uuid: scan}) for s in expand_series(series, {scan.uuid: scan})]

    npt.assert_array_equal([y for _, y, _ in resolved], [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])


def test_resolve_series_of_an_unexpanded_wildcard_draws_its_first_channel():
    """A fit holds the series as typed; it has to resolve to one curve, and the first is it."""
    scan = make_multichannel_scan()
    series = make_series(y_name="detector*")

    _, y, _ = resolve_series(series, {scan.uuid: scan})

    npt.assert_array_equal(y, [2.0, 4.0])
