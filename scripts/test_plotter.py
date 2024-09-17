import matplotlib.pylab as plt

from tavi.data.tavi import TAVI
from tavi.instrument.resolution.cooper_nathans import CN
from tavi.plotter import Plot1DManager, Plot2DManager

# from tavi.instrument.instrument_params.python_dicts.cg4c import cg4c_config_params
# from tests.test_data_folder.test_samples.python_samples.nitio3 import nitio3


def test_plot_scan(tavi):
    datasets = list(tavi.data.keys())[0]
    s1 = tavi.data[datasets]["scan0042"]
    # x, y, xerr, yerr, xlabel, ylabel, title, label
    curve1 = s1.generate_curve(norm_channel="mcu", norm_val=30, rebin_type="grid", rebin_step=0.25)

    s2 = tavi.data[datasets]["scan0043"]
    curve2 = s2.generate_curve(norm_channel="mcu", norm_val=30, rebin_type="grid", rebin_step=0.25)

    p1 = Plot1DManager()
    p1.plot_curve(*curve1)
    p1.plot_curve(*curve2)

    plt.show()


def test_2d_plot(tavi, tas):
    datasets = list(tavi.data.keys())[0]

    scan_list = [tavi.data[datasets][f"scan{i:04}"] for i in range(42, 49, 1)] + [
        tavi.data[datasets][f"scan{i:04}"] for i in range(70, 76, 1)
    ]

    # overplot contour
    sg1 = tavi.generate_scan_group(signals=scan_list, signal_axes=("qh", "en", "detector"))
    contour1 = sg1.generate_contour(rebin_steps=(0.025, 0.1))

    plt2d = Plot2DManager()
    plt2d.zlim = [0, 80]
    plt2d.plot_contour(*contour1)

    # -----------------------------------------------
    # generate 2D rez plot
    # -----------------------------------------------

    projection = ((1, 1, 0), (-1, 1, 0), (0, 0, 1))
    q1 = [-0.5, 0.1, 0.05]  # start, stop, step
    q2 = 0
    q3 = 3
    en = [0, 4, 0.4]  # start, stop, step

    ef = 4.8
    R0 = False

    plt2d.rez_plot(tas, projection, q1, q2, q3, en, ef, R0)


if __name__ == "__main__":
    tavi = TAVI()

    tavi_file_name = "./test_data/tavi_test_exp424.h5"
    tavi.open_file(tavi_file_name)

    tas = CN()
    # tas.load_instrument_from_dicts(cg4c_config_params)
    # tas.load_sam ple(nitio3)

    instrument_config_json_path = "./src/tavi/instrument/instrument_params/cg4c.json"
    sample_json_path = "./test_data/test_samples/nitio3.json"

    tas.load_instrument_params_from_json(instrument_config_json_path)
    tas.load_sample_from_json(sample_json_path)

    test_2d_plot(tavi, tas)

    # test_plot_scan(tavi)

    plt.show()
