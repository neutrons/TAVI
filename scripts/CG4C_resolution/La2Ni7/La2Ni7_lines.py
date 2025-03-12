import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt

from tavi.data.fit import Fit1D
from tavi.data.scan import Scan
from tavi.data.tavi import TAVI
from tavi.instrument.resolution.cooper_nathans import CooperNathans
from tavi.plotter import Plot1D
from tavi.sample import Sample


def plot_s1_scan(hkl, scan_num, fit_range=None):
    s1 = Scan.from_spice(path_to_spice_folder, scan_num=scan_num)
    data = s1.get_data(axes=("del_q, s1", "detector"), norm_to=(1, "mcu"))
    # perform fit
    fit = Fit1D(data, fit_range=fit_range)
    fit.add_signal(model="Gaussian")
    fit.add_background(model="Constant")
    pars = fit.guess()
    # pars_s1["b1_c"].set(min=0)
    result = fit.fit(pars, USE_ERRORBAR=False)
    # plot
    p = Plot1D()
    # data
    p.add_scan(data, fmt="o", label="#{} ({},{},{}) s1 scan".format(scan_num, *hkl))
    # fits
    fwhm = result.params["s1_fwhm"]
    p.add_fit(fit, x=fit.x_to_plot(), label=f"FWHM={fwhm.value:.4f}+/-{fwhm.stderr:.4f}")
    return p


def plot_th2th_scan(hkl, scan_num, fit_range=None):
    th2th = Scan.from_spice(path_to_spice_folder, scan_num=scan_num)
    data = th2th.get_data(axes=("del_q", "detector"), norm_to=(1, "mcu"))
    # perform fit
    fit = Fit1D(data, fit_range=fit_range)
    fit.add_signal(model="Gaussian")
    fit.add_background(model="Constant")
    pars = fit.guess()
    # pars_s1["b1_c"].set(min=0)
    result = fit.fit(pars, USE_ERRORBAR=False)
    # plot
    p = Plot1D()
    # data
    p.add_scan(data, fmt="o", label="#{} ({},{},{}) th2th scan".format(scan_num, *hkl))
    # fits
    fwhm = result.params["s1_fwhm"]
    p.add_fit(fit, x=fit.x_to_plot(), label=f"FWHM={fwhm.value:.4f}+/-{fwhm.stderr:.4f}")
    return p


def plot_tilt_scan(hkl, scan_num, tilt_axis="sgl", fit_range=None):
    tilt = Scan.from_spice(path_to_spice_folder, scan_num=scan_num)
    data = tilt.get_data(axes=("del_q," + tilt_axis, "detector"), norm_to=(1, "mcu"))
    # perform fit
    fit = Fit1D(data, fit_range=fit_range)
    fit.add_signal(model="Gaussian")
    fit.add_background(model="Constant")
    pars = fit.guess()
    pars["b1_c"].set(min=0)
    pars["s1_center"].set(value=0)
    result = fit.fit(pars, USE_ERRORBAR=False)
    # plot
    p = Plot1D()
    # data
    p.add_scan(data, fmt="o", label="#{} ({},{},{}) tilt scan".format(scan_num, *hkl))
    # fits
    fwhm = result.params["s1_fwhm"]
    if fwhm.stderr is None:
        label = f"FWHM={fwhm.value:.4f}"
    else:
        label = f"FWHM={fwhm.value:.4f}+/-{fwhm.stderr:.4f}"
    p.add_fit(fit, x=fit.x_to_plot(), label=label)
    return p


def plot_eng_scan(hkl, scan_num, fit_range=None, fit_incoherent_only=False):
    eng = Scan.from_spice(path_to_spice_folder, scan_num=scan_num)
    data = eng.get_data(axes=("en", "detector"), norm_to=(1, "mcu"))
    # perform fit
    fit = Fit1D(data, fit_range=fit_range)
    fit.add_signal(model="Gaussian")

    if fit_incoherent_only:
        fit.add_background(model="Gaussian")
        pars = fit.guess()
        pars["b1_fwhm"].set(value=0.2)
        pars["b1_amplitude"].set(min=0.0)
        pars["b1_center"].set(value=0, vary=False)
    else:
        fit.add_background(model="Constant")
        pars = fit.guess()
        pars["b1_c"].set(min=0)

    # pars["s1_center"].set(value=0)
    result = fit.fit(pars, USE_ERRORBAR=False)
    # plot
    p = Plot1D()
    # data
    p.add_scan(data, fmt="o", label="#{} ({},{},{}) energy scan".format(scan_num, *hkl))
    # fits
    fwhm = result.params["s1_fwhm"]
    p.add_fit(fit, x=fit.x_to_plot(), label=f"FWHM={fwhm.value:.4f}+/-{fwhm.stderr:.4f}")
    return p


def plot_peaks():

    figs = []
    hkl = (0, 0, 6)
    fig, (ax0, ax1, ax2, ax3) = plt.subplots(ncols=4, sharey=True, figsize=(16, 5))
    plot_s1_scan(hkl, scan_num=16).plot(ax0)
    plot_th2th_scan(hkl, scan_num=17).plot(ax1)
    plot_tilt_scan(hkl, scan_num=18).plot(ax2)
    plot_eng_scan(hkl, scan_num=19).plot(ax3)
    figs.append(fig)

    hkl = (0, 0, 4)
    fig, (ax0, ax1, ax2, ax3) = plt.subplots(ncols=4, sharey=True, figsize=(16, 5))
    plot_s1_scan(hkl, scan_num=20).plot(ax0)
    plot_th2th_scan(hkl, scan_num=21).plot(ax1)
    plot_tilt_scan(hkl, scan_num=22).plot(ax2)
    plot_eng_scan(hkl, scan_num=24).plot(ax3)
    figs.append(fig)

    hkl = (0, 0, 2)
    fig, (ax0, ax1, ax2, ax3) = plt.subplots(ncols=4, sharey=True, figsize=(16, 5))
    plot_s1_scan(hkl, scan_num=25).plot(ax0)
    plot_th2th_scan(hkl, scan_num=26).plot(ax1)
    plot_tilt_scan(hkl, scan_num=27).plot(ax2)
    plot_eng_scan(hkl, scan_num=28).plot(ax3)
    figs.append(fig)

    hkl = (0, 0, 5)
    fig, (ax0, ax1, ax2, ax3) = plt.subplots(ncols=4, sharey=True, figsize=(16, 5))
    plot_s1_scan(hkl, scan_num=29).plot(ax0)
    plot_th2th_scan(hkl, scan_num=30).plot(ax1)
    plot_tilt_scan(hkl, scan_num=31).plot(ax2)
    plot_eng_scan(hkl, scan_num=32, fit_incoherent_only=True).plot(ax3)
    figs.append(fig)

    hkl = (0, 0, 3)
    fig, (ax0, ax1, ax2, ax3) = plt.subplots(ncols=4, sharey=True, figsize=(16, 5))
    plot_s1_scan(hkl, scan_num=33).plot(ax0)
    plot_th2th_scan(hkl, scan_num=34).plot(ax1)
    plot_tilt_scan(hkl, scan_num=35).plot(ax2)
    plot_eng_scan(hkl, scan_num=36, fit_incoherent_only=True).plot(ax3)
    figs.append(fig)

    # (001) scans not measured, (0,0,2) measured again by mistake
    # hkl = (0, 0, 1)
    # fig, (ax0, ax1, ax2) = plt.subplots(ncols=3, sharey=True, figsize=(12, 5))
    # plot_s1_scan(hkl, scan_num=37).plot(ax0)
    # plot_th2th_scan(hkl, scan_num=38).plot(ax1)
    # plot_tilt_scan(hkl, scan_num=39).plot(ax2)
    # figs.append(fig)

    hkl = (1, 0, 0)
    fig, (ax0, ax1, ax2, ax3) = plt.subplots(ncols=4, sharey=True, figsize=(16, 5))
    plot_s1_scan(hkl, scan_num=40).plot(ax0)
    plot_th2th_scan(hkl, scan_num=41).plot(ax1)
    plot_tilt_scan(hkl, tilt_axis="sgu", scan_num=42).plot(ax2)
    # (100) energy scan missing
    plot_eng_scan(hkl, scan_num=146).plot(ax3)
    figs.append(fig)

    hkl = (1, 0, 7)
    fig, (ax0, ax1, ax2) = plt.subplots(ncols=3, sharey=True, figsize=(12, 5))
    plot_s1_scan(hkl, scan_num=44).plot(ax0)
    plot_th2th_scan(hkl, scan_num=45).plot(ax1)
    plot_eng_scan(hkl, scan_num=46).plot(ax2)
    figs.append(fig)

    hkl = (1, 0, 6)
    fig, (ax0, ax1, ax2) = plt.subplots(ncols=3, sharey=True, figsize=(12, 5))
    plot_s1_scan(hkl, scan_num=47).plot(ax0)
    plot_th2th_scan(hkl, scan_num=48).plot(ax1)
    plot_eng_scan(hkl, scan_num=49).plot(ax2)
    figs.append(fig)

    hkl = (1, 0, 5)
    fig, (ax0, ax1, ax2) = plt.subplots(ncols=3, sharey=True, figsize=(12, 5))
    plot_s1_scan(hkl, scan_num=50).plot(ax0)
    plot_th2th_scan(hkl, scan_num=51).plot(ax1)
    plot_eng_scan(hkl, scan_num=52).plot(ax2)
    figs.append(fig)

    hkl = (1, 0, 4)
    fig, (ax0, ax1, ax2) = plt.subplots(ncols=3, sharey=True, figsize=(12, 5))
    plot_s1_scan(hkl, scan_num=53).plot(ax0)
    plot_th2th_scan(hkl, scan_num=54).plot(ax1)
    plot_eng_scan(hkl, scan_num=55).plot(ax2)
    figs.append(fig)

    hkl = (1, 0, 3)
    fig, (ax0, ax1, ax2) = plt.subplots(ncols=3, sharey=True, figsize=(12, 5))
    plot_s1_scan(hkl, scan_num=56).plot(ax0)
    plot_th2th_scan(hkl, scan_num=57).plot(ax1)
    plot_eng_scan(hkl, scan_num=58).plot(ax2)
    figs.append(fig)

    hkl = (1, 0, 2)
    fig, (ax0, ax1, ax2) = plt.subplots(ncols=3, sharey=True, figsize=(12, 5))
    plot_s1_scan(hkl, scan_num=59).plot(ax0)
    plot_th2th_scan(hkl, scan_num=60).plot(ax1)
    plot_eng_scan(hkl, scan_num=61).plot(ax2)
    figs.append(fig)

    hkl = (1, 0, 1)
    fig, (ax0, ax1, ax2) = plt.subplots(ncols=3, sharey=True, figsize=(12, 5))
    plot_s1_scan(hkl, scan_num=62).plot(ax0)
    plot_th2th_scan(hkl, scan_num=63).plot(ax1)
    plot_eng_scan(hkl, scan_num=64).plot(ax2)
    figs.append(fig)

    pdf = matplotlib.backends.backend_pdf.PdfPages("./test_data/CTAX_rez/La2Ni7_CG4C_lines.pdf")
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

    plot_peaks()

    # plt.show()
