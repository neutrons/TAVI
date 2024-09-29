from tavi.data.nxdict import NexusDict
from tavi.data.spice_reader import _create_spicelogs


def test_set_dataset():
    nxsource = NexusDict()
    nxsource.set_dataset(key="name", dataset="HFIR", type="NX_CHAR", EX_required="true")

    pass


def test_spice_scan_to_nxdict():
    path_to_spice_data = "./test_data/exp424/Datafiles/CG4C_exp0424_scan0034.dat"
    spicelogs = _create_spicelogs(path_to_spice_data)
    nxdict = NexusDict(daslogs_dict=spicelogs)

    nxsource = NexusDict()
    nxsource.set_attrs(NX_class="NXsource", EX_required="true")
