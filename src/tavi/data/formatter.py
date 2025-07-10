from typing import Dict

import numpy as np
import scipp as sc
import xarray as xr


def dict_to_sc_array(metadata: Dict, data: Dict):
    """Format SPICE data to scipp.DataArray

    Note:
        Return a DataArray for a single scan with one detector channel,
        return a list of DataArray for HB1 if multiple detector channel exsit in one scan"""

    detecor = data.pop("detector")
    coords = {k: sc.array(dims=["Pt."], values=v) for k, v in data.items()}

    da = sc.DataArray(
        data=sc.array(dims=["Pt."], values=detecor, variances=np.sqrt(detecor)),
        coords=coords,
    )

    return da


def dict_to_xr_data_array(metadata: Dict, data: Dict):
    """Format SPICE data to xarray.DataArray

    Note:
        Return a DataArray for a single scan with one detector channel,
        return a list of DataArray for HB1 if multiple detector channel exsit in one scan"""

    detecor = data.pop("detector")
    coords = {
        k: (
            "Pt.",  # dims
            v,  # data
            {"name": k},  # attrs
        )
        for k, v in data.items()
    }

    da = xr.DataArray(
        data=detecor,
        dims=("Pt.",),
        coords=coords,
        attrs={
            "NX_class": "NXentry",
            "scan": metadata["scan"],
        },
        name="detector",  # lable for y-axis
    )

    return da


def dict_to_xr_dataset(metadata: Dict, data: Dict):
    """Format SPICE data to xarray.DataArray

    Note:
        Return a DataArray for a single scan with one detector channel,
        return a list of DataArray for HB1 if multiple detector channel exsit in one scan"""

    detecor = data.pop("detector")
    coords = {k: ("Pt.", v) for k, v in data.items()}

    ds = xr.Dataset(
        data_vars={
            "signal": (  # tuple of form (dim, data, attrs)
                ["Pt."],  # dim
                [int(d) for d in detecor],  # data
                {"name": "detector", "unit": "counts"},  # attrs
            ),
            "error": (["Pt."], np.sqrt(detecor)),
        },
        coords=coords,
        attrs={
            "scan": metadata["scan"],
        },
    )

    return ds
