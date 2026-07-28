import numpy as np
import pytest

from tavi.library.data.plot import Plot
from tavi.library.data.scan import UUID


def make_plot(**kwargs) -> Plot:
    defaults = dict(
        x=np.array([1.0, 2.0, 3.0]),
        y=np.array([4.0, 5.0, 6.0]),
        err=np.array([0.1, 0.2, 0.3]),
        scan_name="test_scan",
        normalized_by="monitor",
        x_name="qh",
        y_name="en",
        error_name="error",
        source_scan_uuid=UUID(value="scan-001"),
    )
    defaults.update(kwargs)
    return Plot(**defaults)


def test_plot_can_be_created():
    plot = make_plot()

    np.testing.assert_array_equal(plot.x, [1.0, 2.0, 3.0])
    np.testing.assert_array_equal(plot.y, [4.0, 5.0, 6.0])
    np.testing.assert_array_equal(plot.err, [0.1, 0.2, 0.3])
    assert plot.scan_name == "test_scan"
    assert plot.normalized_by == "monitor"
    assert plot.x_name == "qh"
    assert plot.y_name == "en"
    assert plot.error_name == "error"


def test_plot_uuid_auto_generated():
    plot = make_plot()

    assert isinstance(plot.uuid, UUID)
    assert plot.uuid.value != ""


def test_plot_uuid_unique_per_instance():
    p1 = make_plot()
    p2 = make_plot()

    assert p1.uuid.value != p2.uuid.value


def test_plot_normalized_by_can_be_none():
    plot = make_plot(normalized_by=None)

    assert plot.normalized_by is None


def test_plot_stores_numpy_arrays():
    x = np.array([10.0, 20.0])
    y = np.array([30.0, 40.0])
    err = np.zeros(2)

    plot = make_plot(x=x, y=y, err=err)

    assert isinstance(plot.x, np.ndarray)
    assert isinstance(plot.y, np.ndarray)
    assert isinstance(plot.err, np.ndarray)
