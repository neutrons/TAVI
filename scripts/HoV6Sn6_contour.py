import matplotlib.pyplot as plt
from tavi.data.tavi import TAVI

tavi = TAVI()

tavi_file_name = "./test_data/tavi_test_exp1031.h5"
tavi.new_tavi_file(tavi_file_name)

nexus_file_name = "./test_data/nexus_exp1031.h5"
tavi.load_nexus_data_from_disk(nexus_file_name)


dataset = tavi.data["IPTS32912_HB1A_exp1031"]

print(len(dataset))
# 001 s1 scans
scan_list = [dataset[f"scan{i:04}"] for i in range(90, 121, 5)]

sg1 = tavi.generate_scan_group(signals=scan_list, signal_axes=("persistent_field", "s1", "detector"))
contour1 = sg1.generate_contour(rebin_steps=(None, None))
sg1.plot_contour(contour1, cmap="turbo")


plt.show()
