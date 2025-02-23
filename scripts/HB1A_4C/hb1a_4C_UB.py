from tavi.instrument.tas import TAS
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
        angles_list.append(MotorAngles(two_theta=two_theta, omega=omega, chi=chi, phi=phi, sgl=None, sgu=None))
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
        angles_list.append(MotorAngles(two_theta=two_theta, omega=omega, chi=chi, phi=phi, sgl=None, sgu=None))
        h = int(float(line_parts[0]))
        k = int(float(line_parts[1]))
        l = int(float(line_parts[2]))
        hkl_list.append((h, k, l))
    return hkl_list, angles_list


def rucl3():
    sample_json_path = "test_data/IPTS33477_HB1A_exp1012/RuCl3/RuCl3.json"
    sample = Sample.from_json(sample_json_path)
    ub_json = sample.ub_conf.ub_mat
    tas.mount_sample(sample)

    peak1 = Peak(hkl=hkl_list[0], angles=angles_list[0])
    peak2 = Peak(hkl=hkl_list[1], angles=angles_list[1])
    tas.calculate_ub_matrix(peaks=(peak1, peak2))
    print(tas.sample.ub_conf.ub_mat)
    print(peak1)
    angles = tas.calculate_motor_angles(hkl=hkl_list[0])
    print(angles)
    print(peak2)
    angles = tas.calculate_motor_angles(hkl=hkl_list[1])
    print(angles)


def la2ni7():
    sample_json_path = "test_data/IPTS33477_HB1A_exp1012/La2Ni7.json"
    sample = Sample.from_json(sample_json_path)
    ub_json = sample.ub_conf.ub_mat
    tas.mount_sample(sample)

    a1 = MotorAngles(two_theta=-64.7086, omega=-32.3539, chi=26.6593, phi=159.038)
    peak1 = Peak(hkl=(-1, -2, -1), angles=a1)
    a2 = MotorAngles(two_theta=-84.7436, omega=-42.3787, chi=44.5, phi=164.4)
    peak2 = Peak(hkl=(-1, -2, -2), angles=a2)

    ubconf = tas.calculate_ub_matrix(peaks=(peak1, peak2))
    print(f"json UB={ub_json}")
    print(f"Calculated UB={ubconf.ub_mat}")

    a1_cal = tas.calculate_motor_angles(hkl=(-1, -2, -1))
    print(f"Exp agnles for (-1,-2,-1)={a1}")
    print(f"Calculated agnles for (-1,-2,-1)={a1_cal}")

    a2_cal = tas.calculate_motor_angles(hkl=(-1, -2, -2))
    print(f"Exp agnles for (-1,-2,-2)={a2}")
    print(f"Calculated agnles for (-1,-2,-2)={a2_cal}")


if __name__ == "__main__":

    # scan_nums, hkl_list, angles_list = read_macro()
    hkl_list, angles_list = read_reflection_list()

    instrument_config_json_path = "test_data/IPTS33477_HB1A_exp1012/hb1a_4c.json"
    tas = TAS(fixed_ef=14.4643, fixed_ei=14.4503, spice_convention=True)
    tas.load_instrument_params_from_json(instrument_config_json_path)

    la2ni7()

    # def chi_phi(hkl, angles):
    #     q = ub.dot(hkl)
    #     # print(q)
    #     chi = -np.degrees(np.arctan2(q[1], np.sqrt(q[0] ** 2 + q[2] ** 2)))
    #     print(chi)
    #     phi = np.degrees(np.arctan2(q[2], q[0]))
    #     print(phi)

    #     print(f"chi={angles.chi}, phi={angles.phi}.")

    # for i in range(len(hkl_list)):
    #     chi_phi(hkl_list[i], angles_list[i])

    #     hkl = tas.calcvulate_hkl_from_angles(angles_list[i])
    #     print(hkl, hkl_list[i])

    # print(angles_list[0])
    # angles_1 = tas.calculate_motor_angles(hkl=hkl_list[0])
    # print(angles_list[0])
    # print(peak1.angles)
