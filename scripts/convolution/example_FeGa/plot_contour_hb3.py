import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axisartist import Axes

from tavi.data.scan import Scan
from tavi.data.scan_group import ScanGroup
from tavi.instrument.tas import TAS
from tavi.plotter import Plot2D
from tavi.sample import Sample

if __name__ == "__main__":
    # setup HB3
    instrument_config_json_path = "./src/tavi/instrument/instrument_params/hb3.json"
    hb3 = TAS(fixed_ef=14.7)
    hb3.load_instrument_params_from_json(instrument_config_json_path)

    # 28.8% Ga in Fe, phonon T2 mode along (1,1,0), scan #19 - 27
    path_to_spice_folder = "test_data/FeGa_data/exp171"
    fe_ga_29 = Sample.from_scan(Scan.from_spice(path_to_spice_folder, scan_num=27))
    hb3.mount_sample(fe_ga_29)

    scan_list = range(19, 27 + 1)
    scans = [Scan.from_spice(path_to_spice_folder, scan_num=num) for num in scan_list]
    sg = ScanGroup(scans)
    # method 1
    # scan_data_2d = sg.combine_data(axes= ( "qk","en", "detector"),
    #     norm_to= (1, "mcu"),
    #     grid= ((1.1, 1.55, 0.05),(0,12,0.4)))
    # method 2
    scan_data_2d = sg.combine_data_hkle(
        axes=((1, 1, 0), (-1, 1, 0), (0, 0, 1), "en", "detector"),
        norm_to=(1, "mcu"),
        grid=((0.99, 1.01), (0, 0.6, 0.05), (-0.01, 0.01), (0, 15, 0.4)),
    )

    # resolution
    axes = ((1, 1, 0), (-1, 1, 0), (0, 0, 1), "en")
    grid = (1, (0, 0.6, 0.1), 0, (0, 15, 3))

    hkle_list = hb3.generate_hkle(grid, axes)
    rez_list = hb3.cooper_nathans(hkle=hkle_list, axes=axes)
    axes_ellip = [i for i, val in enumerate(grid) if isinstance(val, tuple) and len(val) == 3]

    # plotting
    p = Plot2D()
    p.add_contour(
        scan_data_2d,
        cmap="turbo",
        # norm= Normalize(vmin=0,vmax=10),
        norm=LogNorm(vmin=4e0, vmax=1e3),
    )

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
