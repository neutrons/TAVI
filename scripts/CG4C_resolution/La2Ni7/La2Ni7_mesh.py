import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt

from tavi.data.tavi import TAVI
from tavi.instrument.resolution.cooper_nathans_bak import CooperNathans
from tavi.plotter import Plot2D
from tavi.sample import Sample


def h_en_mesh(scans, name, grid, vmax):
    sg = tavi.group_scans(scans, name=name)
    data = sg.combine_data(
        axes=("qh", "en", "detector"),
        norm_to=(1, "mcu"),
        grid=grid,
    )
    p = Plot2D()
    p.add_contour(data, cmap="turbo", vmin=0, vmax=vmax)
    p.title = sg.name
    return p


def l_en_mesh(scans, name, grid, vmax):
    sg = tavi.group_scans(scans, name=name)
    data = sg.combine_data(
        axes=("ql", "en", "detector"),
        norm_to=(1, "mcu"),
        grid=grid,
    )
    p = Plot2D()
    p.add_contour(data, cmap="turbo", vmin=0, vmax=vmax)
    p.title = sg.name
    return p


def h_l_mesh(scans, name, grid, vmax):
    sg = tavi.group_scans(scans, name=name)
    data = sg.combine_data(
        axes=("qh", "ql", "detector"),
        norm_to=(1, "mcu"),
        grid=grid,
    )
    p = Plot2D()
    p.add_contour(data, cmap="turbo", vmin=0, vmax=vmax)
    p.title = sg.name
    return p


def plot_mesh():
    figs = []
    # -----------------------------------------------------
    # -----------------------------------------------------
    ipts_str = "IPTS34324_CG4C_exp0442"

    # ---------006 -------------

    scan_nums = list(range(87, 107 + 1))
    h0l_scans_006 = [(ipts_str, num) for num in scan_nums]
    p = h_l_mesh(h0l_scans_006, "(0,0,6)", grid=((-0.02, 0.02, 0.002), (5.6, 6.4, 0.01)), vmax=3e3)
    fig, ax = plt.subplots()
    p.plot(ax)
    figs.append(fig)

    scan_nums = list(range(108, 128 + 1))
    h0l_scans_004 = [(ipts_str, num) for num in scan_nums]
    p = h_l_mesh(h0l_scans_004, "(0,0,4)", grid=((-0.02, 0.02, 0.002), (3.6, 4.4, 0.01)), vmax=1e3)
    fig, ax = plt.subplots()
    p.plot(ax)
    figs.append(fig)

    scan_nums = list(range(129, 149 + 1))
    h0l_scans_002 = [(ipts_str, num) for num in scan_nums]
    p = h_l_mesh(h0l_scans_002, "(0,0,2)", grid=((-0.02, 0.02, 0.002), (1.6, 2.4, 0.01)), vmax=1e3)
    fig, ax = plt.subplots()
    p.plot(ax)
    figs.append(fig)

    scan_nums = list(range(151, 201 + 1))
    h0l_scans_100 = [(ipts_str, num) for num in scan_nums]
    p = h_l_mesh(h0l_scans_100, "(1,0,0)", grid=((0.96, 1.04, 0.002), (-0.2, 0.2, 0.02)), vmax=3e3)
    fig, ax = plt.subplots()
    p.plot(ax)
    figs.append(fig)

    # -----------------------------------------------------
    # -----------------------------------------------------

    ipts_str = "IPTS35030_CG4C_exp0445"
    # ------------ 107-------------
    scan_nums = list(range(65, 81 + 1))
    h_scans_107 = [(ipts_str, num) for num in scan_nums]
    p = h_en_mesh(h_scans_107, "(1,0,7)", grid=((0.8, 1.2, 0.002), (-0.5, 0.5, 0.05)), vmax=2e3)
    fig, ax = plt.subplots()
    p.plot(ax)
    figs.append(fig)

    scan_nums = list(range(82, 102 + 1))
    l_scans_h07 = [(ipts_str, num) for num in scan_nums]
    p = l_en_mesh(l_scans_h07, "(1,0,7)", grid=((6.5, 7.5, 0.01), (-0.5, 0.5, 0.05)), vmax=2e3)
    fig, ax = plt.subplots()
    p.plot(ax)
    figs.append(fig)

    # ------------ 004-------------
    scan_nums = list(range(103, 123 + 1))
    h_scans_004 = [(ipts_str, num) for num in scan_nums]
    p = h_en_mesh(h_scans_004, "(0,0,4)", grid=((-0.2, 0.2, 0.004), (-0.5, 0.5, 0.05)), vmax=1e3)
    fig, ax = plt.subplots()
    p.plot(ax)
    figs.append(fig)

    scan_nums = list(range(124, 144 + 1))
    l_scans_004 = [(ipts_str, num) for num in scan_nums]
    p = l_en_mesh(l_scans_004, "(0,0,4)", grid=((3.5, 4.5, 0.01), (-0.5, 0.5, 0.05)), vmax=1e3)
    fig, ax = plt.subplots()
    p.plot(ax)
    figs.append(fig)

    # ------------ 107 hori foc-------------
    scan_nums = list(range(196, 212 + 1))
    h_scans_107 = [(ipts_str, num) for num in scan_nums]
    p = h_en_mesh(h_scans_107, "(1,0,7) w/ horizontal focusing", grid=((0.8, 1.2, 0.002), (-0.5, 0.5, 0.05)), vmax=2e3)
    fig, ax = plt.subplots()
    p.plot(ax)
    figs.append(fig)

    scan_nums = list(range(213, 233 + 1))
    l_scans_h07 = [(ipts_str, num) for num in scan_nums]
    p = l_en_mesh(l_scans_h07, "(1,0,7) w/ horizontal focusing", grid=((6.5, 7.5, 0.01), (-0.5, 0.5, 0.05)), vmax=2e3)
    fig, ax = plt.subplots()
    p.plot(ax)
    figs.append(fig)

    scan_nums = list(range(317, 363 + 1)) + list(range(365, 386 + 1))
    h0l_scans_107 = [(ipts_str, num) for num in scan_nums]
    p = h_l_mesh(h0l_scans_107, "(1,0,7) w/ horizontal focusing", grid=((0.9, 1.1, 0.002), (6.6, 7.5, 0.02)), vmax=2e3)
    fig, ax = plt.subplots()
    p.plot(ax)
    figs.append(fig)

    # ------------ 004 hori foc-------------
    scan_nums = list(range(234, 254 + 1))
    h_scans_004 = [(ipts_str, num) for num in scan_nums]
    p = h_en_mesh(h_scans_004, "(0,0,4) w/ horizontal focusing", grid=((-0.2, 0.2, 0.004), (-0.5, 0.5, 0.05)), vmax=1e3)
    fig, ax = plt.subplots()
    p.plot(ax)
    figs.append(fig)

    scan_nums = list(range(255, 275 + 1))
    l_scans_004 = [(ipts_str, num) for num in scan_nums]
    p = l_en_mesh(l_scans_004, "(0,0,4) w/ horizontal focusing", grid=((3.5, 4.5, 0.01), (-0.5, 0.5, 0.05)), vmax=1e3)
    fig, ax = plt.subplots()
    p.plot(ax)
    figs.append(fig)

    scan_nums = list(range(276, 316 + 1))
    h0l_scans_004 = [(ipts_str, num) for num in scan_nums]
    p = h_l_mesh(
        h0l_scans_004, "(0,0,4) w/ horizontal focusing", grid=((-0.04, 0.04, 0.002), (3.6, 4.4, 0.01)), vmax=1e3
    )
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
    tavi.load_spice_data_from_disk("test_data/CTAX_rez/IPTS-35030/exp445/")
    tavi.load_spice_data_from_disk("test_data/CTAX_rez/IPTS-34324/exp442/")

    plot_mesh()

    # plt.show()
