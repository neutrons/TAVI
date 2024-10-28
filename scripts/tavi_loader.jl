using HDF5

filename = "/Users/x4l/Documents/GitHub/TAVI/test_data/tavi_rez.h5"
fid = h5open(filename, "r")
pdata = fid["processed_data"]
scan0042 = pdata["scan0042"]
rez_mat = read(scan0042, "rez_mat")
