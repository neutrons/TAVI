import matplotlib.pyplot as plt
import numpy as np
import scipp as sc
import xarray as xr
from pytest import fixture

from tavi.data.formatter import dict_to_sc_array, dict_to_xr_data_array, dict_to_xr_dataset
from tavi.data.spice_reader import read_spice_datafile


@fixture
def spice_data():
    path_to_spice_data = "./test_data/exp424/Datafiles/CG4C_exp0424_scan0063.dat"
    metadata, data, *_ = read_spice_datafile(path_to_spice_data)
    path_to_spice_data = "./test_data/exp424/Datafiles/CG4C_exp0424_scan0064.dat"
    metadata2, data2, *_ = read_spice_datafile(path_to_spice_data)
    return ((metadata, data), (metadata2, data2))


def test_sc_array_quick_plot(spice_data):
    sc_array = dict_to_sc_array(*spice_data[0])
    sc_array2 = dict_to_sc_array(*spice_data[1])
    assert isinstance(sc_array, sc.DataArray)

    _, axes = plt.subplots(ncols=2)
    sc_array.plot(
        coords="e",
        ax=axes[0],
    )
    sc_array2.plot(coords="e", ax=axes[0])
    sc_array2.plot(coords="e", ax=axes[1])
    axes[0].legend()
    axes[0].grid(alpha=0.6)
    axes[1].grid(alpha=0.6)
    plt.tight_layout()

    plt.show()


def test_sc_array_plot(spice_data):
    sc_array = dict_to_sc_array(*spice_data[0])

    _, ax = plt.subplots()
    ax.errorbar(
        x=sc_array.coords["e"],
        y=sc_array.data.values,
        yerr=sc_array.data.variances,
        fmt="o",
    )

    plt.show()


def test_xr_data_array_quick_plot(spice_data):
    xr_array = dict_to_xr_data_array(*spice_data[0])
    xr_array2 = dict_to_xr_data_array(*spice_data[1])

    fig, axes = plt.subplots(ncols=2)
    xr_array.plot(ax=axes[0], x="e", marker="o", label=f"{xr_array.attrs['scan']}")
    xr_array2.plot(ax=axes[0], x="e", marker="s", label=f"{xr_array2.attrs['scan']}")
    xr_array.plot(ax=axes[1], x="ei", color="k")
    axes[0].legend()
    axes[0].grid(alpha=0.6)
    axes[1].grid(alpha=0.6)
    plt.tight_layout()
    plt.show()


def test_xr_data_array_plot(spice_data):
    xr_array = dict_to_xr_data_array(*spice_data[0])

    (fig, ax0) = plt.subplots()
    x = xr_array.coords["e"]
    xlab = x.attrs["name"]
    y = xr_array
    yerr = np.sqrt(y)
    ylab = y.name
    ax0.errorbar(x, y, yerr=yerr, fmt="o", label=xr_array.attrs["scan"])
    ax0.set_xlabel(xlab)
    ax0.set_ylabel(ylab)
    ax0.legend()
    ax0.grid(alpha=0.6)
    plt.tight_layout()
    plt.show()


def test_dataset_plotting(spice_data):
    ds = dict_to_xr_dataset(*spice_data[0])
    ds2 = dict_to_xr_dataset(*spice_data[1])
    assert isinstance(ds, xr.Dataset)

    fig, ax = plt.subplots()
    ax.errorbar(
        ds["e"],
        ds["signal"],
        yerr=ds["error"],
        fmt="o",
        label=f"{ds.attrs['scan']}",
    )
    ax.errorbar(
        ds2["e"],
        ds2["signal"],
        yerr=ds2["error"],
        fmt="o",
        label=f"{ds2.attrs['scan']}",
    )
    ax.legend()
    ax.grid(alpha=0.6)
    plt.tight_layout()
    plt.show()


def test_dataset_to_h5(spice_data):
    ds = dict_to_xr_dataset(*spice_data[0])

    root = xr.DataTree(name="root")
    root["ds"] = xr.DataTree(ds, name="ds")
    root.to_netcdf("saved_on_disk.h5", engine="h5netcdf")

    # Read it back
    root = xr.open_datatree("saved_on_disk.h5", engine="h5netcdf")
    print(root)


def test_dataset_renorm(spice_data):
    ds = dict_to_xr_dataset(*spice_data[0])

    signal = ds["signal"]
    norm = ds.coords["mcu"]
    signal_renormed = signal / norm

    _, ax = plt.subplots()
    ax.plot(ds.coords["e"], signal)
    ax.plot(ds.coords["e"], signal_renormed * 30.0)
    plt.show()


def test_sc_xarray_conversion(spice_data):
    xr_array = dict_to_xr_data_array(*spice_data[0])
    assert xr_array.attrs["NX_class"] == "NXentry"
    sc_da = sc.compat.xarray_compat.from_xarray(xr_array)
    assert isinstance(sc_da, sc.DataArray)


def test_xr_data_array_rebin(spice_data):
    da = dict_to_xr_data_array(*spice_data[0])

    bins = np.arange(0, 5, 0.05) - 0.05 / 2  # [0 1,2,3,4]
    da_binned = da.groupby_bins("e", bins=bins, right=False).mean()

    # Extract bin centers from Interval objects
    intervals = da_binned.coords["e_bins"].values  # pandas.Interval objects
    bin_centers = np.array([interval.left + (interval.right - interval.left) / 2 for interval in intervals])

    da_binned = da_binned.assign_coords(e_bins=("e_bins", bin_centers))
    da_binned = da_binned.rename({"e_bins": "e"})

    (fig, ax0) = plt.subplots()

    ax0.plot(da.coords["e"], da, "-o", label="original")
    ax0.plot(da_binned.coords["e"], da_binned, "-o", label="binned")
    ax0.legend()
    ax0.grid(alpha=0.6)
    plt.tight_layout()
    plt.show()
