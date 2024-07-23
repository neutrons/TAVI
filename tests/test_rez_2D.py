import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.axisartist import Axes
from tavi.instrument.instrument_params.takin_test import instrument_params
from test_data_folder.test_samples.sample_test import test_xtal
from tavi.instrument.resolution.cooper_nathans import CN


def rez_plot_gen(instrument_params, xtal, projection, q1, q2, q3, en, ef, R0):

    qe_list = np.empty((4,), dtype=object)
    plot_axes = []
    perp_axes = []
    for idx, qe in enumerate((q1, q2, q3, en)):
        if np.size(qe) > 1:
            plot_axes.append(idx)
            qe_list[idx] = np.arange(qe[0], qe[1] + qe[2] / 2, qe[2])
        else:
            perp_axes.append(idx)
            qe_list[idx] = np.array([qe])
    qe_list[3] = qe_list[3] + ef

    tas = CN()
    tas.load_instrument(instrument_params)
    tas.load_sample(xtal)

    fig = plt.figure(figsize=(10, 6))
    has_axes = False
    p1, p2, p3 = projection

    for q1 in qe_list[0]:
        for q2 in qe_list[1]:
            for q3 in qe_list[2]:
                for ei0 in qe_list[3]:

                    h, k, l = tuple(np.array(p1) * q1 + np.array(p2) * q2 + np.array(p3) * q3)

                    rez = tas.cooper_nathans(
                        ei=ei0,
                        ef=ef,
                        hkl=(h, k, l),
                        projection=projection,
                        R0=R0,
                    )
                    if rez.STATUS:
                        elps_qx_en = rez.generate_ellipse(axes=tuple(plot_axes), PROJECTION=False)
                        if not has_axes:
                            ax = fig.add_subplot(111, axes_class=Axes, grid_helper=elps_qx_en.grid_helper)
                            has_axes = True
                        elps_qx_en.generate_plot(ax, c="black", linestyle="solid")
                        elps_proj_qx_en = rez.generate_ellipse(axes=tuple(plot_axes), PROJECTION=True)
                        elps_proj_qx_en.generate_plot(ax, c="black", linestyle="dashed")

    if has_axes:
        projection = projection + ("en",)
        ax.set_title(
            f"{projection[perp_axes[0]]}={qe_list[perp_axes[0]][0]}, "
            + f"{projection[perp_axes[1]]}={qe_list[perp_axes[1]][0]}"
        )


def test_l_vs_en():

    projection = ((0, 0, 1), (0, -1, 0), (2, -1, 0))

    q1 = [0, 2, 0.5]  # start, stop, step
    q2 = 0
    q3 = 0
    en = [0, 13, 2]  # start, stop, step

    ef = 13.5
    R0 = False

    rez_plot_gen(instrument_params, test_xtal, projection, q1, q2, q3, en, ef, R0)


def test_h_vs_k():

    projection = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    q1 = [0, 2.1, 0.5]  # start, stop, step
    q2 = [0, 2.1, 0.5]  # start, stop, step
    q3 = 0
    en = 5
    ef = 13.5
    R0 = False

    rez_plot_gen(instrument_params, test_xtal, projection, q1, q2, q3, en, ef, R0)


if __name__ == "__main__":
    test_l_vs_en()
    test_h_vs_k()

plt.show()
