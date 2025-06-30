# -*- coding: utf-8 -*

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axisartist import Axes

from tavi.data.tavi import TAVI
from tavi.instrument.tas import TAS
from tavi.plotter import Plot2D
from tavi.sample import Sample
from tavi.utilities import labels_from_projection


def test_plot2d():
    instrument_config_json_path = "./test_data/IPTS9879_HB1A_exp978/hb1a.json"
    tas = TAS(fixed_ei=14.450292, convention="Spice")
    tas.load_instrument_params_from_json(instrument_config_json_path)

    tavi = TAVI()
    path_to_spice_folder = "./test_data/IPTS9879_HB1A_exp978/exp978/"
    tavi.load_spice_data_from_disk(path_to_spice_folder)

    # -------------- Si (111) 40'-40'-40'-80' -----------------
    scans = list(range(824, 848 + 1))

    sg1 = tavi.group_scans(scans, name="Si (111) 40'-40'-40'-80'")
    si_111_1 = sg1.combine_data(
        axes=("qh", "ql", "detector"),
        norm_to=(1, "mcu"),
        grid=((0.97, 1.03, 0.0025), (0.97, 1.03, 0.0025)),
    )
    si_111_2 = sg1.combine_data(
        axes=("s1", "s2", "detector"),
        norm_to=(1, "mcu"),
        grid=((-24.5, -21.5, 0.1), (-46.5, -43, 0.1)),
    )

    p1 = Plot2D()
    p1.add_contour(si_111_1, cmap="turbo", vmin=0, vmax=2.5e4)
    p1.title = sg1.name
    p1.ylim = [0.97, 1.03]
    p1.xlim = [0.97, 1.03]

    p2 = Plot2D()
    p2.add_contour(si_111_2, cmap="turbo", vmin=0, vmax=2.5e4)
    p2.title = sg1.name
    # p2.ylim = [0.97, 1.03]
    # p2.xlim = [0.97, 1.03]

    fig, axes = plt.subplots(ncols=2)
    im1 = p1.plot(axes[0])
    im2 = p2.plot(axes[1])
    # fig.colorbar()
    plt.show()


# TODO
def test_plot2d_with_resolution():
    # load data
    tavi = TAVI("./test_data/tavi_exp424.h5")
    scan_list = list(range(42, 49, 1)) + list(range(70, 76, 1))

    sg = tavi.group_scans(scan_list, name="dispH")
    scan_data_2d = sg.combine_data(
        axes=("qh", "en", "detector"),
        norm_to=(1, "mcu"),
        grid=(0.025, (-0.5, 4.5, 0.1)),
    )
    # load experimental parameters
    instrument_config_json_path = "./src/tavi/instrument/instrument_params/cg4c.json"
    tas = TAS(fixed_ef=4.8)
    tas.load_instrument_params_from_json(instrument_config_json_path)

    sample_json_path = "./test_data/test_samples/nitio3.json"
    sample = Sample.from_json(sample_json_path)
    tas.mount_sample(sample)

    # calculate resolution ellipses
    hkle_list = [(qh, qh, 3, en) for qh in np.arange(-0.5, 0.15, 0.05) for en in np.arange(0, 4.1, 0.4)]
    axes = ((1, 1, 0), (0, 0, 1), (1, -1, 0), "en")
    rez_list = tas.cooper_nathans(hkle=hkle_list, axes=axes)

    # generate plot
    p = Plot2D()
    p.add_contour(scan_data_2d, cmap="turbo", vmax=2)

    for rez in rez_list:
        e_co = rez.get_ellipse(axes=(3, 0), PROJECTION=False)
        e_inco = rez.get_ellipse(axes=(3, 0), PROJECTION=True)
        p.add_reso(e_co, c="k", linestyle="solid")
        p.add_reso(e_inco, c="k", linestyle="dashed")

    fig = plt.figure()
    ax = fig.add_subplot(111, axes_class=Axes)

    im = p.plot(ax)
    fig.colorbar(im, ax=ax)
    plt.show()


def test_making_labels_from_projection():
    label = labels_from_projection()
    assert label == ("(H, 0, 0) (r.l.u.)", "(0, K, 0) (r.l.u.)", "(0, 0, L) (r.l.u.)", "E (meV)")

    label = labels_from_projection(axes=None)
    assert label == ("Q_para (A^-1)", "Q_perp (A^-1)", "Q_up (A^-1)", "E (meV)")

    label = labels_from_projection(axes=((1, 1, 0), (0, 0, 1), (1, -1, 0), "en"))
    assert label == ("(H, H, 0) (r.l.u.)", "(0, 0, L) (r.l.u.)", "(K, -K, 0) (r.l.u.)", "E (meV)")

    label = labels_from_projection(axes=((1.0, 1.0, 0.0), (0, 0, 1), "en", (1, -1, 0)))
    assert label == ("(H, H, 0) (r.l.u.)", "(0, 0, L) (r.l.u.)", "E (meV)", "(K, -K, 0) (r.l.u.)")

    label = labels_from_projection(axes=((0, 0, 1), (1, -1, 0), (2, 2, 0), "en"))
    assert label == ("(0, 0, L) (r.l.u.)", "(H, -H, 0) (r.l.u.)", "(2K, 2K, 0) (r.l.u.)", "E (meV)")
