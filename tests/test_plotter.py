# -*- coding: utf-8 -*

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axisartist import Axes

from tavi.data.tavi import TAVI
from tavi.instrument.resolution.cooper_nathans import CN
from tavi.plotter import Plot2D
from tavi.sample.xtal import Xtal


def test_plot2d():

    # load data
    tavi = TAVI("./test_data/tavi_exp424.h5")
    scan_list = list(range(42, 49, 1)) + list(range(70, 76, 1))

    sg = tavi.combine_scans(scan_list, name="dispH")
    scan_data_2d = sg.get_data(
        axes=("qh", "en", "detector"),
        norm_to=(1, "mcu"),
        grid=(0.025, (-0.5, 4.5, 0.1)),
    )
    # load experimental parameters
    instrument_config_json_path = "./src/tavi/instrument/instrument_params/cg4c.json"
    tas = CN(SPICE_CONVENTION=False)
    tas.load_instrument_params_from_json(instrument_config_json_path)

    sample_json_path = "./test_data/test_samples/nitio3.json"
    sample = Xtal.from_json(sample_json_path)
    tas.mount_sample(sample)

    # calculate resolution ellipses
    R0 = False
    hkl_list = [(qh, qh, 3) for qh in np.arange(-0.5, 0.1, 0.05)]
    ef = 4.8
    ei_list = [e + ef for e in np.arange(0, 4.1, 0.4)]
    projection = ((1, 1, 0), (0, 0, 1), (1, -1, 0))
    rez_list = tas.cooper_nathans(hkl_list=hkl_list, ei=ei_list, ef=ef, projection=projection, R0=R0)

    # genreate plot
    p = Plot2D()
    im = p.add_contour(scan_data_2d, cmap="turbo", vmax=1)

    for rez in rez_list:
        e_co = rez.get_ellipse(axes=(0, 3), PROJECTION=False)
        e_inco = rez.get_ellipse(axes=(0, 3), PROJECTION=True)
        p.add_reso(e_co, c="k", linestyle="solid")
        p.add_reso(e_inco, c="k", linestyle="dashed")

    fig = plt.figure()
    ax = fig.add_subplot(111, axes_class=Axes, grid_helper=p.grid_helper)

    im = p.plot(ax)
    fig.colorbar(im, ax=ax)
    plt.show()
