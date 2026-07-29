import numpy.testing as npt

from tavi.library.data.plot import Plot, PlotSeries
from tavi.library.data.plot_resolution import resolve_series, scans_for_plots
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


def test_resolve_series_error_is_sqrt_abs_y_over_two():
    scan = make_scan()
    series = make_series()

    _, y, err = resolve_series(series, {scan.uuid: scan})

    npt.assert_allclose(err, [1.0, 1.5, 2.0])


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


def test_resolve_series_normalizes_err_the_same_way_as_y():
    scan = make_scan(data={"qh": [1.0, 2.0, 3.0], "en": [4.0, 9.0, 16.0], "monitor": [2.0, 3.0, 4.0]})
    series = make_series(normalized_by="monitor", normalized_by_value=1.0)

    _, _, err = resolve_series(series, {scan.uuid: scan})

    npt.assert_allclose(err, [1.0 / 2.0, 1.5 / 3.0, 2.0 / 4.0])


def test_scans_for_plots_returns_only_referenced_scans():
    scan_a = make_scan(uuid_val="scan-a")
    scan_b = make_scan(uuid_val="scan-b")
    plot = Plot(series=[make_series(uuid_val="scan-a")])

    result = scans_for_plots([plot], {scan_a.uuid: scan_a, scan_b.uuid: scan_b})

    assert set(result) == {scan_a.uuid}
