import xarray as xr

# Example coordinates and detector data
# x = [1, 2, 3, 4, 5, 6, 7]
# detector = [10.5, 11.0, 9.8, 12.3, 13.1, 14.0, 15.2]

# da = xr.DataArray(
#     data=detector,
#     dims="points",
#     coords={
#         "points": x,
#         "x": ("points", [x0 * 2 for x0 in x]),
#     },
#     name="detector",
# )

# print(da)


# # Define bin edges
# bins = np.arange(0, 8, 2)  # [0, 2, 4, 6]

# # Group by x bins
# binned = da.groupby_bins("points", bins=bins, right=False).mean()

# print(binned)


# Create datasets
ds_root = xr.Dataset({"a": ("x", [1, 2, 3])})
ds_child = xr.Dataset({"b": ("x", [4, 5, 6])})

# Create DataTree
tree = xr.DataTree(ds_root, name="root")
tree["child"] = xr.DataTree(ds_child, name="child")

# Write to HDF5/NetCDF file
tree.to_netcdf("tree_output.h5", engine="h5netcdf")

# Read it back
tree_loaded = xr.open_datatree("tree_output.h5", engine="h5netcdf")

print(tree_loaded)
