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


def test_scan_group_contour(tavi):

    print(len(tavi.data))
    scan_list = [tavi.data[f"scan{i:04}"] for i in range(42, 49, 1)] + [
        tavi.data[f"scan{i:04}"] for i in range(70, 76, 1)
    ]

    sg1 = tavi.generate_scan_group(signals=scan_list, signal_axes=("qh", "en", "detector"))
    contour1 = sg1.generate_contour(rebin_steps=(0.025, 0.1))
    sg1.plot_contour(contour1, cmap="turbo", vmax=80)

    sg2 = tavi.generate_scan_group(signals=scan_list, signal_axes=("en", "qh", "detector"))
    contour2 = sg2.generate_contour(rebin_steps=(0.1, 0.025))
    sg2.plot_contour(contour2, cmap="turbo", vmax=80)

    sg3 = tavi.generate_scan_group(
        signals=scan_list,
        signal_axes=("qk", "ei", "detector"),
    )
    contour3 = sg3.generate_contour(rebin_steps=(None, 0.1), norm_channel="mcu", norm_val=30)
    sg3.plot_contour(contour3, cmap="turbo", vmax=40)

    plt.show()


def test_scan_group_waterfall(tavi):

    print(len(tavi.data))
    scan_list = [tavi.data[f"scan{i:04}"] for i in range(42, 49, 1)] + [
        tavi.data[f"scan{i:04}"] for i in range(70, 76, 1)
    ]

    sg = tavi.generate_scan_group(signals=scan_list, signal_axes=("en", "detector", "qh"))

    wf1 = sg.generate_waterfall(norm_channel="mcu", norm_val=600)
    sg.plot_waterfall(wf1, ylim=[0, 1.5e4], xlim=[0, 6], fmt="-o", shifts=1e3)

    plt.show()


if __name__ == "__main__":

    tavi = TAVI_Data()

    tavi_file_name = "./tests/test_data_folder/tavi_test.h5"
    tavi.new_tavi_file(tavi_file_name)

    nexus_file_name = "./tests/test_data_folder/nexus_exp424.h5"
    tavi.load_tavi_data_from_disk(nexus_file_name)

    # test_plot_scan(tavi)
    # test_scan_group_contour(tavi)
    test_scan_group_waterfall(tavi)
