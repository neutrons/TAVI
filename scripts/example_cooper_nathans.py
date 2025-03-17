import matplotlib.pylab as plt
import numpy as np

from tavi.instrument.resolution.cooper_nathans_bak import CooperNathans
from tavi.sample.xtal import Xtal

np.set_printoptions(floatmode="fixed", precision=4)


def test_copper_nathans_localQ(tas_params):
    tas, ei, ef, hkl, _, r0 = tas_params

    rez = tas.cooper_nathans(ei=ei, ef=ef, hkl=hkl, projection=None, r0=r0)
    rez.plot()


def test_copper_nathans_hkl(tas_params):
    tas, ei, ef, hkl, _, R0 = tas_params

    rez = tas.cooper_nathans(ei=ei, ef=ef, hkl=hkl, R0=R0)
    rez.plot()


def test_copper_nathans_projection(tas_params):
    tas, ei, ef, hkl, projection, R0 = tas_params

    rez = tas.cooper_nathans(ei=ei, ef=ef, hkl=hkl, projection=projection, R0=R0)

    rez.plot()


def test_cooper_nathans_compare_3():
    tas = CooperNathans()

    instrument_config_json_path = "./src/tavi/instrument/instrument_params/takin_test.json"
    sample_json_path = "./test_data/test_samples/sample_test.json"

    tas.load_instrument_params_from_json(instrument_config_json_path)
    tas.mount_sample(Xtal.from_json(sample_json_path))

    if tas.sample.ub_mat is None:
        peak_list = [(0, 0, 2), (0, 2, 0)]
        angles_list = [
            (-51.530388, -45.220125, -0.000500, -2.501000),
            (-105.358735, 17.790125, -0.000500, -2.501000),
        ]
        tas.find_ub(peaks=peak_list, angles=angles_list, ei=13.500172, ef=13.505137)

    ei = 13.5
    ef = 13.5
    hkl = (0, 0, 2)
    # projection = ((1, 0, 0), (-1, 2, 0), (0, 0, 1))
    projection = ((0, 0, 1), (0, -1, 0), (2, -1, 0))
    r0 = False

    tas_params = (tas, ei, ef, hkl, projection, r0)

    test_copper_nathans_localQ(tas_params)
    test_copper_nathans_hkl(tas_params)
    test_copper_nathans_projection(tas_params)


def test_cooper_nathans_CTAX():
    tas = CooperNathans()

    instrument_config_json_path = "./src/tavi/instrument/instrument_params/cg4c.json"
    sample_json_path = "./test_data/test_samples/nitio3.json"

    tas.load_instrument_params_from_json(instrument_config_json_path)
    tas.mount_sample(sample_json_path)

    # ----------- reset UB -----------------------------------
    # print(f"u={tas.sample.u}")
    # print(f"v={tas.sample.v}")
    # tas.sample.ub_mat = tas.sample.uv_to_ub_matrix(u=[1, 1, 0], v=[0, 0, 1])
    # print(f"u={tas.sample.u}")
    # print(f"v={tas.sample.v}")
    # ----------------------------------------------------------

    ei = 4.8
    ef = 4.8
    hkl = (0, 0, 3)

    projection = ((1, 1, 0), (0, 0, 1), (1, -1, 0))
    r0 = True

    tas_params = (tas, ei, ef, hkl, projection, r0)

    test_copper_nathans_localQ(tas_params)
    test_copper_nathans_hkl(tas_params)
    test_copper_nathans_projection(tas_params)


if __name__ == "__main__":
    test_cooper_nathans_compare_3()

    test_cooper_nathans_CTAX()

    plt.show()
