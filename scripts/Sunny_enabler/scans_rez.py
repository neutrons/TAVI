from tavi.data.nexus_builder import NXdataset, NXentry
from tavi.data.nexus_entry import NexusEntry
from tavi.data.tavi import TAVI
from tavi.instrument.resolution.cooper_nathans_bak import CooperNathans
from tavi.sample.xtal import Xtal

instrument_config_json_path = "./src/tavi/instrument/instrument_params/cg4c.json"
tas = CooperNathans(SPICE_CONVENTION=False)
tas.load_instrument_params_from_json(instrument_config_json_path)

sample_json_path = "./test_data/test_samples/nitio3.json"
sample = Xtal.from_json(sample_json_path)
tas.mount_sample(sample)

tavi = TAVI()
tavi.load_spice_data_from_disk("./test_data/exp424")

scan_list = list(range(42, 49, 1)) + list(range(70, 76, 1))

for scan_num in scan_list:
    scan = tavi.get_scan(scan_num)

    scan_data = scan.get_data(axes=("en", "detector"))

    ei_list = scan.data["ei"]
    ef_list = scan.data["ef"]
    qh_list = scan.data["qh"]
    qk_list = scan.data["qk"]
    ql_list = scan.data["ql"]

    # projection = ((1, 1, 0), (0, 0, 1), (1, -1, 0))
    R0 = True

    rez_mat_list = []
    rez_r0_list = []

    for i in range(len(ei_list)):

        rez = tas.rez(
            ei=ei_list[i],
            ef=ef_list[i],
            hkl_list=(qh_list[i], qk_list[i], ql_list[i]),
            # projection=projection,
            R0=R0,
        )
        rez_mat_list.append(rez.mat)
        rez_r0_list.append(rez.r0)

    rez_entry = NXentry(
        qh=NXdataset(ds=qh_list),
        qk=NXdataset(ds=qk_list),
        ql=NXdataset(ds=ql_list),
        ei=NXdataset(ds=ei_list),
        ef=NXdataset(ds=ef_list),
        detector=NXdataset(ds=scan_data.y),
        err=NXdataset(ds=scan_data.err),
        rez_mat=NXdataset(ds=rez_mat_list),
        rez_r0=NXdataset(ds=rez_r0_list),
        # projection=NXdataset(ds=projection),
    )

    tavi.processed_data.update(NexusEntry._dict_to_nexus_entry({scan.name: rez_entry}))


file_path = "./test_data/tavi_rez_HKLE.h5"
tavi.save(file_path)
