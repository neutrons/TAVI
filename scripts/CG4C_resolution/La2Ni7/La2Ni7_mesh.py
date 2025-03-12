import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt

from tavi.data.tavi import TAVI
from tavi.instrument.resolution.cooper_nathans import CooperNathans
from tavi.plotter import Plot2D
from tavi.sample import Sample


def h_mesh(scans, name, grid, vmax):

    sg = tavi.combine_scans(scans, name=name)
    data = sg.get_data(
        axes=("qh", "en", "detector"),
        norm_to=(1, "mcu"),
        grid=grid,
    )
    p = Plot2D()
    p.add_contour(data, cmap="turbo", vmin=0, vmax=vmax)
    p.title = sg.name
    return p


def l_mesh(scans, name, grid, vmax):

    sg = tavi.combine_scans(scans, name=name)
    data = sg.get_data(
        axes=("ql", "en", "detector"),
        norm_to=(1, "mcu"),
        grid=grid,
    )
    p = Plot2D()
    p.add_contour(data, cmap="turbo", vmin=0, vmax=vmax)
    p.title = sg.name
    return p


def plot_mesh():
    figs = []
    # ------------ 107-------------
    h_scans_107 = list(range(65, 81 + 1))
    p = h_mesh(h_scans_107, "(1,0,7)", grid=((0.8, 1.2, 0.002), (-0.5, 0.5, 0.05)), vmax=2e3)
    fig, ax = plt.subplots()
    p.plot(ax)
    figs.append(fig)

    l_scans_h07 = list(range(82, 102 + 1))
    p = l_mesh(l_scans_h07, "(1,0,7)", grid=((6.5, 7.5, 0.01), (-0.5, 0.5, 0.05)), vmax=2e3)
    fig, ax = plt.subplots()
    p.plot(ax)
    figs.append(fig)

    # ------------ 004-------------
    h_scans_004 = list(range(103, 123 + 1))
    p = h_mesh(h_scans_004, "(0,0,4)", grid=((-0.2, 0.2, 0.004), (-0.5, 0.5, 0.05)), vmax=1e3)
    fig, ax = plt.subplots()
    p.plot(ax)
    figs.append(fig)

    l_scans_004 = list(range(124, 144 + 1))
    p = l_mesh(l_scans_004, "(0,0,4)", grid=((3.5, 4.5, 0.01), (-0.5, 0.5, 0.05)), vmax=1e3)
    fig, ax = plt.subplots()
    p.plot(ax)
    figs.append(fig)

    # ------------ 107 hori foc-------------
    h_scans_107 = list(range(196, 205 + 1))
    p = h_mesh(h_scans_107, "(1,0,7) w/ horizontal focusing", grid=((0.8, 1.2, 0.002), (-0.5, 0.5, 0.05)), vmax=2e3)
    fig, ax = plt.subplots()
    p.plot(ax)
    figs.append(fig)

    pdf = matplotlib.backends.backend_pdf.PdfPages("./test_data/CTAX_rez/La2Ni7_CG4C_QEmesh.pdf")
    for f in figs:
        pdf.savefig(f)
    pdf.close()


if __name__ == "__main__":

    instrument_config_json_path = "./test_data/CTAX_rez/cg4c.json"
    cg4c = CooperNathans(fixed_ef=4.8, spice_convention=True)
    cg4c.load_instrument_params_from_json(instrument_config_json_path)

    sample_json_path = "./test_data/CTAX_rez/La2Ni7.json"
    cg4c.mount_sample(Sample.from_json(sample_json_path))

    # ------------------------ load data ------------------------
    tavi = TAVI()
    path_to_spice_folder = "test_data/CTAX_rez/IPTS-35030/exp445/"
    tavi.load_spice_data_from_disk(path_to_spice_folder)

    plot_mesh()

    # plt.show()
