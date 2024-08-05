import matplotlib.pylab as plt
from tavi.tavi_data.tavi_data import TAVI_Data
from tavi.plotter import Plot2DManager
from tavi.instrument.resolution.cooper_nathans import CN
from tavi.instrument.instrument_params.python_dicts.cg4c import cg4c_config_params
from tests.test_data_folder.test_samples.python_samples.nitio3 import nitio3

if __name__ == "__main__":

    tavi = TAVI_Data()

    tavi_file_name = "./tests/test_data_folder/tavi_test_exp424.h5"
    tavi.open_tavi_file(tavi_file_name)

    dataset = tavi.data["IPTS32124_CG4C_exp424"]

    scan_list = [dataset[f"scan{i:04}"] for i in range(42, 49, 1)] + [dataset[f"scan{i:04}"] for i in range(70, 76, 1)]

    # overplot contour
    sg1 = tavi.generate_scan_group(signals=scan_list, signal_axes=("qh", "en", "detector"))
    contour1 = sg1.generate_contour(rebin_steps=(0.025, 0.1))

    x, y, z, _, _, xlabel, ylabel, zlabel, title = contour1

    plt2d = Plot2DManager()
    # ax = plt2d.fig.add_subplot(111)
    p = plt2d.ax.pcolormesh(x, y, z, shading="auto", cmap="turbo", vmax=80)
    plt2d.fig.colorbar(p, ax=plt2d.ax)
    plt2d.ax.set_title(title)
    plt2d.ax.set_xlabel(xlabel)
    plt2d.ax.set_ylabel(ylabel)
    plt2d.ax.grid(alpha=0.6)

    # if xlim is not None:
    #     ax.set_xlim(left=xlim[0], right=xlim[1])
    # if ylim is not None:
    #     ax.set_ylim(bottom=ylim[0], top=ylim[1])

    # generate 2D rez plot

    tas = CN()
    tas.load_instrument_from_dicts(cg4c_config_params)
    tas.load_sample(nitio3)

    projection = ((1, 1, 0), (-1, 1, 0), (0, 0, 1))
    q1 = [-0.5, 0.1, 0.05]  # start, stop, step
    q2 = 0
    q3 = 3
    en = [0, 4, 0.4]  # start, stop, step

    ef = 4.8
    R0 = False

    plt2d.rez_plot_gen(tas, projection, q1, q2, q3, en, ef, R0)

    plt.show()
