import matplotlib.pyplot as plt

from tavi.data.tavi import TAVI


def test_plot_scan(tavi):
    print(len(tavi.data))
    datasets = list(tavi.data.keys())[0]
    print(datasets)
    s = tavi.data[datasets]["scan0042"]
    s.plot_curve()
    # s.plot_curve(norm_channel="mcu", norm_val=30, rebin_type="grid", rebin_step=0.25)
    s.plot_curve(rebin_type="grid", rebin_step=0.25)
    s.plot_curve(norm_channel="time", norm_val=30, rebin_type="grid", rebin_step=0.25)
    s.plot_curve(rebin_type="tol", rebin_step=0.25)
    s.plot_curve(norm_channel="time", norm_val=30, rebin_type="tol", rebin_step=0.25)

    plt.show()


def test_scan_group_contour_exp710():
    tavi = TAVI()

    tavi_file_name = "./tests/test_data_folder/tavi_test_exp710.h5"
    tavi.new_file(tavi_file_name)

    nexus_file_name = "./tests/test_data_folder/nexus_exp710.h5"
    tavi.get_nexus_data_from_disk(nexus_file_name)

    scan_list = (
        [tavi.data[f"scan{i:04}"] for i in range(214, 225, 1)]
        + [tavi.data[f"scan{i:04}"] for i in range(279, 282, 1)]
        + [tavi.data[f"scan{i:04}"] for i in range(50, 60, 1)]
        # + [tavi.data[f"scan{i:04}"] for i in range(91, 101, 1)]
    )

    sg1 = tavi.generate_scan_group(signals=scan_list, signal_axes=("qh", "en", "detector"))
    # contour1 = sg1.generate_contour(rebin_steps=(0.1, 0.5))
    contour1 = sg1.generate_contour(rebin_steps=(None, None))
    sg1.plot_contour(contour1, cmap="turbo", vmax=40, ylim=[0, 70])

    sg2 = tavi.generate_scan_group(signals=scan_list, signal_axes=("ei", "qk", "detector"))
    contour2 = sg2.generate_contour()
    sg2.plot_contour(contour2, cmap="turbo", vmax=40)

    plt.show()


def test_scan_group_contour(tavi):
    dataset = tavi.data["IPTS32124_CG4C_exp424"]

    print(len(dataset))
    scan_list = [dataset[f"scan{i:04}"] for i in range(42, 49, 1)] + [dataset[f"scan{i:04}"] for i in range(70, 76, 1)]

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
    contour3 = sg3.generate_contour(rebin_steps=(0.025, 0.1), norm_channel="mcu", norm_val=30)
    sg3.plot_contour(contour3, cmap="turbo", vmax=40)

    plt.show()


def test_scan_group_waterfall(tavi):
    print(len(tavi.data))
    scan_list = [tavi.data[f"scan{i:04}"] for i in range(42, 49, 1)] + [
        tavi.data[f"scan{i:04}"] for i in range(70, 76, 1)
    ]

    # sg1 = tavi.generate_scan_group(signals=scan_list, signal_axes=("qh", "en", "detector"))
    # wf1 = sg1.generate_contour(rebin_steps=(0.025, 0.1), norm_channel="mcu", norm_val=120)
    # sg1.plot_waterfall(wf1, ylim=[0, 1.4e3], xlim=[0, 6.5], fmt="o", shifts=1e2)

    sg2 = tavi.generate_scan_group(signals=scan_list, signal_axes=("en", "qh", "detector"))
    wf2 = sg2.generate_contour(rebin_steps=(0.1, 0.025), norm_channel="mcu", norm_val=120)
    sg2.plot_waterfall(wf2, ylim=[0, 1.4e5], xlim=[-0.6, 0.2], fmt="-o", shifts=1e2)

    plt.show()


if __name__ == "__main__":
    tavi = TAVI()

    tavi_file_name = "./test_data/tavi_test_exp424.h5"
    tavi.open_file(tavi_file_name)

    # test_scan_group_contour(tavi)

    # test_scan_group_waterfall(tavi)

    # test_scan_group_contour_exp710()

    test_plot_scan(tavi)
