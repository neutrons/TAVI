from tavi.instrument.tas import TAS
from tavi.sample import Sample
from tavi.ub_algorithm import ub_matrix_to_uv
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
        angles_list.append(MotorAngles(two_theta=two_theta, omega=omega, chi=chi, phi=phi, sgl=None, sgu=None))
    return scan_nums, hkl_list, angles_list


def read_reflection_list():

    filename = "test_data/IPTS33477_HB1A_exp1012/MnWO4/Refections_022425.list"
    with open(filename, encoding="utf-8") as f:
        all_content = f.readlines()

    hkl_list = []
    angles_list = []

    for i in range(2, len(all_content) - 2):
        line_parts = all_content[i].split("\t")

        two_theta = float(line_parts[0])
        omega = float(line_parts[1])
        chi = float(line_parts[2])
        phi = float(line_parts[3])
        angles_list.append(
            MotorAngles(
                two_theta=two_theta,
                omega=omega,
                chi=chi,
                phi=phi,
                sgl=None,
                sgu=None,
            )
        )
        h = int(float(line_parts[4]))
        k = int(float(line_parts[5]))
        l = int(float(line_parts[6]))
        hkl_list.append((h, k, l))
    return hkl_list, angles_list


def mnwo4():
    sample_json_path = "test_data/IPTS33477_HB1A_exp1012/MnWO4/MnWO4.json"
    sample = Sample.from_json(sample_json_path)
    ub_json = sample.ub_conf.ub_mat
    tas.mount_sample(sample)
    u, v = ub_matrix_to_uv(spice_to_mantid(ub_json))
    print(u, v)
    print(f"json UB=\n{ub_json}")
    print(sample)

    peak1 = Peak(hkl=hkl_list[0], angles=angles_list[0])
    peak2 = Peak(hkl=hkl_list[1], angles=angles_list[1])
    peak3 = Peak(hkl=hkl_list[2], angles=angles_list[2])

    ubconf_2 = tas.calculate_ub_matrix(peaks=(peak1, peak2))
    print(f"Calculated from two peaks UB=\n{ubconf_2.ub_mat}")
    print(tas.sample)

    ubconf_3 = tas.calculate_ub_matrix(peaks=(peak1, peak2, peak3))
    print(f"Calculated from three peaks UB=\n{ubconf_3.ub_mat}")
    print(tas.sample)

    # peaks = tuple([Peak(hkl=hkl_list[i], angles=angles_list[i]) for i in range(len(hkl_list))])
    # ubconf_all = tas.calculate_ub_matrix(peaks=peaks)
    # print(f"Calculated from all peaks UB=\n{ubconf_all.ub_mat}")
    # print(tas.sample)

    for i in range(len(hkl_list)):
        angle_cal = tas.calculate_motor_angles(hkl=hkl_list[i])
        if not (angles_list[i] == angle_cal):
            print(f"Experiment agnles for {hkl_list[i]} are {angles_list[i]}")
            print(f"Calculated agnles for {hkl_list[i]} are {angle_cal}")


if __name__ == "__main__":

    # scan_nums, hkl_list, angles_list = read_macro()
    hkl_list, angles_list = read_reflection_list()

    instrument_config_json_path = "test_data/IPTS33477_HB1A_exp1012/hb1a_4c.json"
    tas = TAS(fixed_ef=14.4643)  # , fixed_ei=14.4503, spice_convention=True)
    tas.load_instrument_params_from_json(instrument_config_json_path)

    mnwo4()
