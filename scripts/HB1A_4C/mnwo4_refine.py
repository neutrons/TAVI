import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
import numpy as np

from tavi.data.fit import Fit1D
from tavi.data.scan import Scan
from tavi.data.tavi import TAVI
from tavi.instrument.tas import TAS
from tavi.plotter import Plot1D
from tavi.sample import Sample
from tavi.utilities import MotorAngles, Peak


def read_macro():
    scan_macro = "test_data/IPTS33477_HB1A_exp1012/exp1012/Macros/MnWO4_RT_overnight2.macro"
    with open(scan_macro, encoding="utf-8") as f:
        all_content = f.readlines()

    hkl_list = []
    angles_list = []
    scan_nums = []
    for i in range(int(len(all_content) / 3)):
        scan_nums.append(i + 662)
        line1 = all_content[3 * i].split(" ")
        qh = int(float(line1[-4][1:]))
        qk = int(float(line1[-3]))
        ql = int(float(line1[-2]))
        hkl_list.append((qh, qk, ql))
        line2 = all_content[3 * i + 1].split(" ")
        two_theta = float(line2[2])
        omega = float(line2[4])
        chi = float(line2[6])
        phi = float(line2[8])
        angles_list.append(MotorAngles(two_theta=two_theta, omega=omega, chi=chi, phi=phi, sgl=None, sgu=None))
    return scan_nums, hkl_list, angles_list


def check_ub():
    """
    -1.3309484e-01   1.2999069e-01  -3.6163741e-02
    1.5660985e-01   1.1374652e-01   6.3277766e-03
    2.5325472e-02  -2.0245888e-02  -1.9664176e-01
    """

    peak1 = Peak(hkl=hkl_list[0], angles=angles_list[0])
    peak2 = Peak(hkl=hkl_list[1], angles=angles_list[1])
    peak3 = Peak(hkl=hkl_list[3], angles=angles_list[3])

    ubconf_2 = hb1a_4c.calculate_ub_matrix(peaks=(peak1, peak2))
    print(f"Calculated from two peaks UB=\n{ubconf_2.ub_mat}")
    print(hb1a_4c.sample)

    ubconf_3 = hb1a_4c.calculate_ub_matrix(peaks=(peak1, peak2, peak3))
    print(f"Calculated from three peaks UB=\n{ubconf_3.ub_mat}")
    print(hb1a_4c.sample)

    peaks = tuple([Peak(hkl=hkl_list[i], angles=angles_list[i]) for i in range(len(hkl_list))])
    ubconf_all = hb1a_4c.calculate_ub_matrix(peaks=peaks)
    print(f"Calculated from all peaks UB=\n{ubconf_all.ub_mat}")
    print(hb1a_4c.sample)

    for i in range(len(hkl_list)):
        angle_cal = hb1a_4c.calculate_motor_angles(hkl=hkl_list[i])
        if not (angles_list[i] == angle_cal):
            print(f"Experiment agnles for {hkl_list[i]} are {angles_list[i]}")
            print(f"Calculated agnles for {hkl_list[i]} are {angle_cal}")


def analyze_s1_scan_in_q(hkl, s1_scan, fit_range=None):
    s1 = Scan.from_spice(path_to_spice_folder, scan_num=s1_scan)

    # --------------- analyze in del_q ----------------
    scan_s1 = s1.get_data(axes=("del_q", "detector"), norm_to=(1, "mcu"))
    # perform fit
    scan_s1_fit = Fit1D(scan_s1, fit_range)
    scan_s1_fit.add_signal(model="Gaussian")
    scan_s1_fit.add_background(model="Constant")
    pars_s1 = scan_s1_fit.guess()
    # pars_s1["b1_c"].set(min=0)
    result_s1 = scan_s1_fit.fit(pars_s1, USE_ERRORBAR=False)
    # print(scan_s1_fit.result.fit_report())
    # print(f"Fit {scan2}")

    p = Plot1D()
    # data
    p.add_scan(scan_s1, fmt="o", label="#{} ({},{},{}) s1 scan".format(s1_scan, *hkl))
    # fits
    fwhm = result_s1.params["s1_fwhm"]
    p.add_fit(
        scan_s1_fit,
        x=scan_s1_fit.x_to_plot(),
        label=f"FWHM={fwhm.value:.4f}+/-{fwhm.stderr:.4f}",
    )

    p.ylim = (-np.max(scan_s1.y) * 0.1, np.max(scan_s1.y) * 1.3)

    return (result_s1, p, np.mean(s1.data.get("q")))


def analyze_s1_scan_in_omega(hkl, s1_scan, fit_range=None):
    s1 = Scan.from_spice(path_to_spice_folder, scan_num=s1_scan)

    # --------------- analyze in del_q ----------------
    scan_s1 = s1.get_data(axes=("s1", "detector"), norm_to=(1, "mcu"))
    # perform fit
    scan_s1_fit = Fit1D(scan_s1, fit_range)
    scan_s1_fit.add_signal(model="Gaussian")
    scan_s1_fit.add_background(model="Constant")
    pars_s1 = scan_s1_fit.guess()
    # pars_s1["b1_c"].set(min=0)
    result_s1 = scan_s1_fit.fit(pars_s1, USE_ERRORBAR=False)
    # print(scan_s1_fit.result.fit_report())
    # print(f"Fit {scan2}")

    p = Plot1D()
    # data
    p.add_scan(scan_s1, fmt="o", label="#{} ({},{},{}) s1 scan".format(s1_scan, *hkl))
    # fits
    fwhm = result_s1.params["s1_fwhm"]
    p.add_fit(
        scan_s1_fit,
        x=scan_s1_fit.x_to_plot(),
        label=f"FWHM={fwhm.value:.4f}+/-{fwhm.stderr:.4f}",
    )

    p.ylim = (-np.max(scan_s1.y) * 0.1, np.max(scan_s1.y) * 1.3)

    return (result_s1, p, np.mean(s1.data.get("q")))


def plot_s1_nuclear_peaks():
    two_theta_list = []
    q_list = []
    rez_list = []
    exp_s1_omega = []
    exp_s1_q = []
    for i, hkl in enumerate(hkl_list):
        s1_q, p1, q = analyze_s1_scan_in_q(hkl, scan_nums[i])

        rez = hb1a_4c.cooper_nathans(hkl=hkl, axes=None)
        p1.add_reso_bar(
            pos=s1_q,
            fwhm=rez.coh_fwhms(axis=1),
            c="C3",
            label=f"Resolution FWHM={rez.coh_fwhms(axis=1):.04f}",
        )

        s1_omega, p2, _ = analyze_s1_scan_in_omega(hkl, scan_nums[i])

        q_list.append(q)
        rez_list.append(rez)
        exp_s1_omega.append(s1_omega)
        exp_s1_q.append(s1_q)
        # two_theta_list.append(two_theta)

        # make plot
        fig, (ax0, ax1) = plt.subplots(ncols=2, sharey=True, figsize=(10, 5))
        p1.plot(ax0)
        p2.plot(ax1)
        # plt.show()

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()

    return (hkl_list, rez_list, exp_s1_omega, exp_s1_q)


def export_intensity(file_name):
    # 690 (-1,1,-1) and 726 (1,-1,1) not good
    lorentz_factor_list = [
        rez.r0 * np.sqrt((rez.mat[0, 0] * rez.mat[1, 1] - rez.mat[0, 1] * rez.mat[1, 0]) / rez.mat[1, 1] / (2 * np.pi))
        for rez in analysis[1]
    ]
    intensity_list = [
        exp_s1_q.params["s1_amplitude"].value / lorentz_factor_list[i] for i, exp_s1_q in enumerate(analysis[3])
    ]
    err_list = [
        exp_s1_q.params["s1_amplitude"].stderr / lorentz_factor_list[i] for i, exp_s1_q in enumerate(analysis[3])
    ]

    with open(file_name, "w") as f:
        f.write("Single crystal data of MnWO4 (hb1a)\n")
        f.write("(3i5,2f8.2,i4,3f8.2)\n")
        f.write("2.38000  0   0\n")

        for i, hkl in enumerate(hkl_list):
            if hkl == (-1, 1, -1) or hkl == (1, -1, 1):
                continue
            h, k, l = hkl
            f.write(f"{h:5d}{k:5d}{l:5d}{intensity_list[i]:8.2f}{err_list[i]:8.2f}   1\n")

    f.close()


if __name__ == "__main__":
    scan_nums, hkl_list, angles_list = read_macro()

    instrument_config_json_path = "test_data/IPTS33477_HB1A_exp1012/hb1a_4c.json"
    # using Spice convention by default
    hb1a_4c = TAS(fixed_ef=14.4643)
    hb1a_4c.load_instrument_params_from_json(instrument_config_json_path)

    sample_json_path = "test_data/IPTS33477_HB1A_exp1012/MnWO4/MnWO4.json"
    mnwo4 = Sample.from_json(sample_json_path)
    print("initial lattice parameters")
    print(mnwo4)

    hb1a_4c.mount_sample(mnwo4)
    peaks = tuple([Peak(hkl=hkl_list[i], angles=angles_list[i]) for i in range(len(hkl_list))])
    ubconf_all = hb1a_4c.calculate_ub_matrix(peaks=peaks)

    print("lattice from refined UB matrix:")
    print(hb1a_4c.sample)

    # check_ub()
    tavi = TAVI()
    path_to_spice_folder = "test_data/IPTS33477_HB1A_exp1012/exp1012/"
    tavi.load_spice_data_from_disk(path_to_spice_folder)

    pdf = matplotlib.backends.backend_pdf.PdfPages("./test_data/IPTS33477_HB1A_exp1012/MnWO4_nuclear_peaks.pdf")
    analysis = plot_s1_nuclear_peaks()
    pdf.close()

    int_file = "./test_data/IPTS33477_HB1A_exp1012/MnWO4_lorentz_removed.int"
    export_intensity(int_file)
