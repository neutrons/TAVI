from concurrent.futures import ProcessPoolExecutor
from functools import partial
from time import time

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axisartist import Axes
from matplotlib.colors import Normalize,LogNorm

from tavi.data.scan import Scan
from tavi.instrument.resolution.convolution import convolution
from tavi.instrument.tas import TAS
from tavi.plotter import Plot2D
from tavi.sample import Sample
from tavi.data.scan_group import ScanGroup

if __name__ == "__main__":
    # -------------------------------------------------------
    # setup HB1A
    instrument_config_json_path = "./src/tavi/instrument/instrument_params/hb1a.json"
    hb1a = TAS(fixed_ei=14.632980)
    hb1a.load_instrument_params_from_json(instrument_config_json_path)
   
    # 10.8% Ga in Fe, phonon T2 mode along (1,1,0), scan #27 - 34 
    path_to_spice_folder = "test_data/FeGa_data/exp266"
    fe_ga_11 = Sample.from_scan(Scan.from_spice(path_to_spice_folder, scan_num=34))
    hb1a.mount_sample(fe_ga_11)

    scan_list = range(27, 34 + 1)

    scans = [Scan.from_spice(path_to_spice_folder, scan_num=num) for num in scan_list]
    sg = ScanGroup(scans)
    # method 1
    # scan_data_2d = sg.combine_data(axes= ( "qk","en", "detector"),
    #     norm_to= (1, "mcu"),
    #     grid= ((1.1, 1.55, 0.05),(-23,0,0.6)))
    scan_data_2d = sg.combine_data_hkle(
        axes = ((1, 1, 0), (-1, 1, 0), (0, 0, 1), "en","detector"),
        norm_to= (1, "mcu"),
        grid= ( (0.99,1.01), (0, 0.6, 0.05),(-0.01,0.01), (-23, 0, 0.6)),)

     # resolution
    axes = ((1, 1, 0), (-1, 1, 0), (0, 0, 1), "en")
    grid = (1,(0,0.6,0.1), 0, (-23, 0, 3))

    hkle_list = hb1a.generate_hkle(grid, axes)
    rez_list = hb1a.cooper_nathans(hkle=hkle_list, axes=axes)
    axes_ellip = [i for i, val in enumerate(grid) if isinstance(val, tuple) and len(val) == 3]
    
    p = Plot2D()
    p.add_contour(scan_data_2d,cmap="turbo",
        #norm= Normalize(vmin=0,vmax=10),
        norm= LogNorm(vmin=1e-1, vmax=1e1))
    

    for rez in filter(None, rez_list):
        e_co = rez.get_ellipse(axes=axes_ellip, PROJECTION=False)
        e_inco = rez.get_ellipse(axes=axes_ellip, PROJECTION=True)
        p.add_reso(e_co, c="k", linestyle="solid")
        p.add_reso(e_inco, c="k", linestyle="dashed")
    
    fig = plt.figure()
    ax = fig.add_subplot(111, axes_class=Axes)
    ax.set_title(p.title)
    im = p.plot(ax)
    fig.colorbar(im, ax=ax)
   
    plt.show()