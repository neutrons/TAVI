import matplotlib.pyplot as plt
import numpy as np

from tavi.data.convfit import ConvFit1D
from tavi.data.fit import Fit1D
from tavi.data.scan import Scan
from tavi.instrument.tas import TAS
from tavi.plotter import Plot1D
from tavi.sample import Sample
from tavi.utilities import MotorAngles, Peak


def gaussian(x, cen, sigma, amp, bkgd=0):
    prefactor = amp / np.sqrt(2 * np.pi) / sigma
    y = np.exp(-((x - cen) ** 2) / 2 / sigma**2) * prefactor
    return y + bkgd


if __name__ == "__main__":
    # load instrument parameters
    instrument_config_json_path = "test_data/IPTS9879_HB1A_exp978/hb1a_La2Ni7.json"
    ei = 14.450292
    ef = 14.443601
    tas = TAS(fixed_ef=ef, convention="Spice")
    tas.load_instrument_params_from_json(instrument_config_json_path)
    # load sample parameters
    sample_json_path = "test_data/IPTS9879_HB1A_exp978/La2Ni7.json"
    sample = Sample.from_json(sample_json_path)
    ub_json = sample.ub_conf.ub_mat
    tas.mount_sample(sample)

    # -------------- check UB calculation -----------------

    angles1 = MotorAngles(two_theta=-101.396853, omega=-48.004475, sgl=-0.770162, sgu=1.477665)
    peak1 = Peak(hkl=(0, 0, 16), angles=angles1)
    angles2 = MotorAngles(two_theta=-56.150124, omega=64.624337, sgl=-0.770162, sgu=1.477665)
    peak2 = Peak(hkl=(1, 1, 0), angles=angles2)

    tas.calculate_ub_matrix(peaks=(peak1, peak2))
    assert np.allclose(tas.sample.ub_conf.ub_mat, ub_json, atol=1e-4)
    # load scan data
    path_to_spice_folder = "test_data/IPTS9879_HB1A_exp978/exp978/"
    scan_num = 152  # th2th
    hkl = (0, 0, 12)
    th2th = Scan.from_spice(path_to_spice_folder, scan_num=scan_num)
    scan_th2th = th2th.get_data(axes=("q", "detector"), norm_to=(1, "mcu"))

    # perform fit
    scan_th2th_fit = Fit1D(scan_th2th)
    scan_th2th_fit.add_signal(model="Gaussian")
    scan_th2th_fit.add_background(model="Constant")
    pars_th2th = scan_th2th_fit.guess()
    pars_th2th["b1_c"].set(min=0)
    result_th2th = scan_th2th_fit.fit(pars_th2th, USE_ERRORBAR=False)
    print(scan_th2th_fit.result.fit_report())

    tas.monochromator.mosaic_h = 30
    tas.analyzer.mosaic_h = 30
    tas.collimators.h_pre_mono = 50
    rez = tas.cooper_nathans(hkl=hkl, en=ei - ef, axes=None)

    # resolution fwhm
    x_th2th = scan_th2th_fit.result.params["s1_center"].value
    components_th2th = result_th2th.eval_components(result_th2th.params, x=x_th2th)
    y_th2th = components_th2th["s1_"] / 2 + components_th2th["b1_"]
    # make plot
    p1 = Plot1D()
    # data
    p1.add_scan(scan_th2th, fmt="o", label="#{} ({},{},{}) th2th scan".format(scan_num, *hkl))
    p1.add_reso_bar(
        pos=(x_th2th, y_th2th),
        fwhm=rez.coh_fwhms(axis=0),
        c="k",
        label=f"Resolution FWHM={rez.coh_fwhms(axis=0):.04f}",
    )
    p1.add_fit(
        scan_th2th_fit,
        x=scan_th2th_fit.x_to_plot(),
        label=f"Gaussian FWHM={result_th2th.params['s1_fwhm'].value:.4f}+/-{result_th2th.params['s1_fwhm'].stderr:.4f}",
    )

    # perform fit convoluted with resolution
    rez_params = th2th.get_coherent_fwhm(tas=tas, projection=None, axis=0)
    scan_th2th_fit_rez = ConvFit1D(scan_th2th, rez_params)
    scan_th2th_fit_rez.add_signal(model="Gaussian")
    # scan_th2th_fit_rez.add_background(model="Constant")
    pars_th2th_rez = scan_th2th_fit_rez.guess()
    # pars_th2th_rez["s1_sigma"].value = 0.01
    # pars_th2th_rez["s1_amplitude"].value = 1e6
    result_th2th_rez = scan_th2th_fit_rez.fit(pars_th2th_rez, USE_ERRORBAR=False)
    print(result_th2th_rez.fit_report())

    # fits
    fwhm = result_th2th_rez.params["s1_fwhm"]
    p1.add_fit(
        scan_th2th_fit_rez,
        label="Fit Convoluted w/ Resolution",
    )
    p1.add_fit(
        scan_th2th_fit_rez,
        DECONV=True,
        x=scan_th2th_fit_rez.x_to_plot(num_of_pts=500),
        label=f"Intrinsic FWHM={fwhm.value:.4f}+/-{fwhm.stderr:.4f}",
    )

    # plotting
    fig, ax = plt.subplots(figsize=(8, 6))
    p1.plot(ax)
    plt.show()
