import matplotlib.pyplot as plt
import numpy as np

from tavi.data.fit import Fit1D
from tavi.data.scan import Scan
from tavi.instrument.resolution.cooper_nathans import CooperNathans
from tavi.plotter import Plot1D
from tavi.sample import Sample
from tavi.utilities import MotorAngles, Peak, spice_to_mantid


def read_macro():
    scan_macro = "test_data/IPTS33477_HB1A_exp1012/exp1012/Macros/MnWO4_RT_overnight2.macro"
    with open(scan_macro, encoding="utf-8") as f:
        all_content = f.readlines()

    hkl_list = []
    angles_list = []
    scan_nums = []
    for i in range(int(len(all_content) / 3)):
        scan_nums.append(3 * i + 662)
        line1 = all_content[3 * i].split(" ")
        h = int(float(line1[-4][1:]))
        k = int(float(line1[-3]))
        l = int(float(line1[-2]))
        hkl_list.append((h, k, l))
        line2 = all_content[3 * i + 1].split(" ")
        two_theta = float(line2[2])
        omega = float(line2[4])
        chi = float(line2[6])
        phi = float(line2[8])
        angles_list.append(
            MotorAngles(
                two_theta=two_theta,
                omega=omega,
                chi=chi,
                phi=phi,
                sgl=0,
                sgu=0,
            )
        )
    return scan_nums, hkl_list, angles_list


def read_reflection_list():
    # filename = "test_data/IPTS33477_HB1A_exp1012/Reflections_012825.list"
    filename = "test_data/IPTS33477_HB1A_exp1012/RuCl3/fit_observ.dat"
    with open(filename, encoding="utf-8") as f:
        all_content = f.readlines()

    hkl_list = []
    angles_list = []

    for i in range(2, len(all_content) - 2):
        line_parts = all_content[i].split(" ")
        two_theta = float(line_parts[3])
        omega = float(line_parts[4])
        chi = float(line_parts[5])
        phi = float(line_parts[6][:-1])
        angles_list.append(
            MotorAngles(
                two_theta=two_theta,
                omega=omega,
                chi=chi,
                phi=phi,
                sgl=0,
                sgu=0,
            )
        )
        h = int(float(line_parts[0]))
        k = int(float(line_parts[1]))
        l = int(float(line_parts[2]))
        hkl_list.append((h, k, l))
    return hkl_list, angles_list


def analyze_in_angles(hkl, scans, fit_ranges):
    scan1, scan2 = scans
    fit_range1, fit_range2 = fit_ranges

    # ------------ resolution -------------
    rez = tas.rez(hkl_list=hkl, ei=ei, ef=ef, R0=False, projection=None)
    # ------------------------- s1 -------------------------
    # print(scan2)
    s1 = Scan.from_spice(path_to_spice_folder, scan_num=scan2)
    scan_s1 = s1.get_data(axes=("s1", "detector"), norm_to=(1, "mcu"))
    # perform fit
    scan_s1_fit = Fit1D(scan_s1, fit_range2)
    scan_s1_fit.add_signal(model="Gaussian")
    scan_s1_fit.add_background(model="Constant")
    pars_s1 = scan_s1_fit.guess()
    # pars_s1["b1_c"].set(min=0)
    result_s1 = scan_s1_fit.fit(pars_s1, USE_ERRORBAR=False)
    # print(scan_s1_fit.result.fit_report())

    p2 = Plot1D()
    # data
    p2.add_scan(scan_s1, fmt="o", label="#{} ({},{},{}) s1 scan".format(scan2, *hkl))
    # fits
    if result_s1.success:
        y = result_s1.params["s1_fwhm"].value
        if (err := result_s1.params["s1_fwhm"].stderr) is None:
            err = 0
        p2.add_fit(
            scan_s1_fit,
            x=scan_s1_fit.x_to_plot(),
            label=f"FWHM={y:.4f}+/-{err:.4f}",
        )
    # resolution
    x_s1 = result_s1.params["s1_center"].value
    components_s1 = result_s1.eval_components(result_s1.params, x=x_s1)
    y_s1 = components_s1["s1_"] / 2 + components_s1["b1_"]
    # convert q to s1
    q_avg = np.mean(s1.data.get("q"))
    fwhm_s1 = np.rad2deg(rez.coh_fwhms(axis=1) / q_avg)
    p2.add_reso_bar(
        pos=(x_s1, y_s1),
        fwhm=fwhm_s1,
        c="C3",
        label=f"Resolution FWHM={fwhm_s1:.04f}",
    )

    p2.ylim = (-np.max(scan_s1.y) * 0.1, np.max(scan_s1.y) * 1.3)

    # make plot
    fig, axes = plt.subplots(ncols=2, sharey=True, figsize=(10, 5))
    p2.plot(axes[1])

    return (np.mean(s1.data.get("q")), result_s1, rez, theta)


if __name__ == "__main__":

    # scan_nums, hkl_list, angles_list = read_macro()
    hkl_list, angles_list = read_reflection_list()

    instrument_config_json_path = "test_data/IPTS33477_HB1A_exp1012/hb1a_4c.json"
    tas = CooperNathans(fixed_ef=14.4643, fixed_ei=14.4503)
    tas.load_instrument_params_from_json(instrument_config_json_path)

    sample_json_path = "test_data/IPTS33477_HB1A_exp1012/RuCl3/RuCl3.json"
    sample = Sample.from_json(sample_json_path)
    ub_json = sample.ub_conf.ub_mat
    tas.mount_sample(sample)

    peak1 = Peak(hkl=hkl_list[0], angles=angles_list[0])
    peak2 = Peak(hkl=hkl_list[1], angles=angles_list[1])
    tas.calculate_ub_matrix(peaks=(peak1, peak2))
    print(tas.sample.ub_conf.ub_mat)
    ub = spice_to_mantid(tas.sample.ub_conf.ub_mat)
    print(ub)

    def chi_phi(hkl, angles):
        q = ub.dot(hkl)
        # print(q)
        chi = -np.degrees(np.arctan2(q[1], np.sqrt(q[0] ** 2 + q[2] ** 2)))
        print(chi)
        phi = np.degrees(np.arctan2(q[2], q[0]))
        print(phi)

        print(f"chi={angles.chi}, phi={angles.phi}.")

    for i in range(len(hkl_list)):
        chi_phi(hkl_list[i], angles_list[i])

        hkl = tas.calcvulate_hkl_from_angles(angles_list[i])
        print(hkl, hkl_list[i])

    # print(angles_list[0])
    # angles_1 = tas.calculate_motor_angles(hkl=hkl_list[0])
    # print(angles_list[0])
    # print(peak1.angles)
