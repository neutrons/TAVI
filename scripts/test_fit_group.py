import matplotlib.pyplot as plt

from tavi.data.plotter import Plot1DManager, Plot2DManager
from tavi.data.tavi import TAVI


def test_fit_group(tavi):
    scan_list1 = [dataset[f"scan{i:04}"] for i in range(90, 121, 5)] + [
        dataset[f"scan{i:04}"] for i in range(127, 128, 5)
    ]
    sg1 = tavi.generate_scan_group(signals=scan_list1, signal_axes=("persistent_field", "s1", "detector"))

    # --------------- plot contour ---------------
    contour1 = sg1.generate_contour(rebin_steps=(None, None), norm_channel="mcu", norm_val=1)
    plt2d1 = Plot2DManager()
    plt2d1.zlim = [0, 5e2]
    plt2d1.title = "(0,0,1)"
    plt2d1.plot_contour(*contour1)

    # x, y, xerr, yerr, xlabel, ylabel, title, label = curve1
    # f1 = Fit(x=x, y=y, fit_range=(0.0, 4))
    # f1.add_background(values=(0.7,))
    # f1.add_signal(
    #     values=(None, 3.5, None),
    #     vary=(True, True, True),
    # )
    # f1.add_signal(
    #     model="Gaussian",
    #     values=(None, 0, None),
    #     vary=(True, False, True),
    #     mins=(0, 0, 0.1),
    #     maxs=(None, 0.1, 0.3),
    # )
    # f1.perform_fit()

    p1 = Plot1DManager()
    # p1.plot_curve(*curve1)
    # p1.plot_curve(f1.x_plot, f1.y_plot, fmt="-")


if __name__ == "__main__":
    tavi = TAVI()

    tavi_file_name = "./test_data/tavi_test_exp1031.h5"
    tavi.new(tavi_file_name)

    nexus_file_name = "./test_data/nexus_exp1031.h5"
    tavi.load_nexus_data_from_disk(nexus_file_name)
    dataset = tavi.data["IPTS32912_HB1A_exp1031"]

    test_fit_group(tavi)

    plt.show()
