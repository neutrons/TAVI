import scipp as sc
import scippnexus as snx

# Create example DataArray
data = sc.array(dims=["x"], values=[1, 2, 3], unit="m")
da = sc.DataArray(data)

# # Save
# sc.io.hdf5.save_hdf5(da, "test_data.h5")

# # Load
# da_loaded = sc.io.hdf5.load_hdf5("test_data.h5")
# print(da_loaded)


filename = "./test_data/IPTS32124_CG4C_exp0424/scan0002.h5"
ins = snx.load(filename, root="scan0002/instrument")
with snx.File(filename) as f:
    scan = f["scan0002"]

pass
