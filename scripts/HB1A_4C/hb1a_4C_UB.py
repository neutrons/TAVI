from tavi.instrument.tas import TAS
from tavi.sample import Sample
from tavi.utilities import MotorAngles, Peak


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
        angles_list.append(MotorAngles(two_theta=two_theta, omega=omega, chi=chi, phi=phi, sgl=None, sgu=None))
        h = int(float(line_parts[4]))
        k = int(float(line_parts[5]))
        l = int(float(line_parts[6]))
        hkl_list.append((h, k, l))
    return hkl_list, angles_list


def mnwo4():
    sample_json_path = "test_data/IPTS33477_HB1A_exp1012/MnWO4/MnWO4.json"
    sample = Sample.from_json(sample_json_path)
    ub_json = sample.ub_conf.ub_mat
    hb1a_4c.mount_sample(sample)
    u, v = hb1a_4c.uv
    print(f"u={u}, v={v}")
    print(f"json UB=\n{ub_json}")
    print(sample)

    peak1 = Peak(hkl=hkl_list[0], angles=angles_list[0])
    peak2 = Peak(hkl=hkl_list[1], angles=angles_list[1])
    peak3 = Peak(hkl=hkl_list[2], angles=angles_list[2])

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

    print(f"Calculated from all peaks in Mantid convention UB=\n{ubconf_all._ub_mat}")

    for i in range(len(hkl_list)):
        angle_cal = hb1a_4c.calculate_motor_angles(hkl=hkl_list[i])
        if not (angles_list[i] == angle_cal):
            print(f"Experiment agnles for {hkl_list[i]} are {angles_list[i]}")
            print(f"Calculated agnles for {hkl_list[i]} are {angle_cal}")


if __name__ == "__main__":
    # scan_nums, hkl_list, angles_list = read_macro()
    hkl_list, angles_list = read_reflection_list()

    instrument_config_json_path = "test_data/IPTS33477_HB1A_exp1012/hb1a_4c.json"
    # using Spice convention by default
    hb1a_4c = TAS(fixed_ef=14.4643)
    hb1a_4c.load_instrument_params_from_json(instrument_config_json_path)

    mnwo4()
