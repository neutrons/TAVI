import matplotlib.pyplot as plt
from tavi.data.fit import Fit
from tavi.data.tavi import TAVI
from tavi.plotter import Plot1DManager


def test_fit_scan(tavi):
    p1 = Plot1DManager()

    datasets = list(tavi.data.keys())[0]
    s1 = tavi.data[datasets]["scan0042"]
    curve1 = s1.generate_curve(norm_channel="mcu", norm_val=30, rebin_type="grid", rebin_step=0.25)

    x, y, xerr, yerr, xlabel, ylabel, title, label = curve1
    f1 = Fit(x=x, y=y, fit_range=(0.0, 4))
    f1.add_background()
    f1.add_signal()
    f1.add_signal(model="Gaussian", values=(0, None, None), vary=(False, True, True))
    f1.perform_fit()

    p1.plot_curve(*curve1)
    p1.plot_curve(f1.x_plot, f1.y_plot, fmt="-")

    # s2 = tavi.data[datasets]["scan0043"]
    # curve2 = s2.generate_curve(norm_channel="mcu", norm_val=30, rebin_type="grid", rebin_step=0.25)
    # p1.plot_curve(*curve2)

    plt.show()


if __name__ == "__main__":
    tavi = TAVI()

    tavi_file_name = "./test_data/tavi_test_exp424.h5"
    tavi.open_tavi_file(tavi_file_name)

    test_fit_scan(tavi)
