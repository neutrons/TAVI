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
    print(peak1)
    print(peak2)
    print(peak3)

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

    check_ub()
