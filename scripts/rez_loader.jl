using HDF5
filename = "./test_data/tavi_rez.h5"
fid = h5open(filename, "r")
pdata = fid["processed_data"]
scan0042 = pdata["scan0042"]
rez_mat = read(scan0042, "rez_mat")

scan0001 = fid["data"]["IPTS32124_CG4C_exp0424"]["scan0001"]
read(scan0001["sample"], "unit_cell")