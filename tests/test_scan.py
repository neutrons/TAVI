import matplotlib.pyplot as plt
from tavi.tavi_data.tavi_data import TAVI_Data


def test_plot_scan(tavi):

    print(len(tavi.data))
    s = tavi.data["scan0025"]
    s.plot_curve()
    # s.plot_curve(norm_channel="mcu", norm_val=30, rebin_type="grid", rebin_step=0.25)
    s.plot_curve(rebin_type="grid", rebin_step=0.25)
    s.plot_curve(norm_channel="time", norm_val=30, rebin_type="grid", rebin_step=0.25)
    s.plot_curve(rebin_type="tol", rebin_step=0.25)
    s.plot_curve(norm_channel="time", norm_val=30, rebin_type="tol", rebin_step=0.25)

    plt.show()


def test_scan_group(tavi):

    print(len(tavi.data))
    scan_list = [tavi.data[f"scan{i:04}"] for i in range(42, 49, 1)]

    sg1 = tavi.generate_scan_group(signals=scan_list, signal_x="qh", signal_y="en", signal_z="detector")
    sg1.plot_image()
    sg2 = tavi.generate_scan_group(signals=scan_list, signal_y="qh", signal_x="en", signal_z="detector")
    sg2.plot_image()

    plt.show()


if __name__ == "__main__":

    tavi = TAVI_Data()

    tavi_file_name = "./tests/test_data_folder/tavi_test.h5"
    tavi.new_tavi_file(tavi_file_name)

    nexus_file_name = "./tests/test_data_folder/nexus_exp424.h5"
    tavi.load_tavi_data_from_disk(nexus_file_name)

    # test_plot_scan(tavi)

    test_scan_group(tavi)
