def spice_to_data_array(file_name: str) -> Union[xr.DataArray, list[xr.DataArray]]:
    """Format SPICE data to xarray.DataArray

    Note:
        Return a DataArray for a single scan with one detector channel,
        return a list of DataArray for HB1 if multiple detector channel exsit in one scan"""

    metadata, data, unrecognized, error_msg = read_spice_datafile(file_name)
    return xr.DataArray()
