from tavi.instrument.resolution.cooper_nathans import CooperNathans
from tavi.sample.xtal import Xtal

instrument_config_json_path = "./src/tavi/instrument/instrument_params/cg4c.json"
tas = CooperNathans(SPICE_CONVENTION=False)
tas.load_instrument_params_from_json(instrument_config_json_path)

sample_json_path = "./test_data/test_samples/nitio3.json"
sample = Xtal.from_json(sample_json_path)
tas.mount_sample(sample)

ei = 4.8
ef = 4.8
hkl = (0, 0, 3)

projection = ((1, 1, 0), (0, 0, 1), (1, -1, 0))
R0 = False

tas_params = (tas, ei, ef, hkl, projection, R0)
rez_list = tas.rez(hkl_list=[(0, 0, 3), (0, 0, -3)], ei=[ei, ei + 1], ef=ef, projection=None)
print(rez_list[0].mat)
