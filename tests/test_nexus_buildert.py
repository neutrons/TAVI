import numpy as np

from tavi.data.nexus_builder import NXdataset, NXentry, spice_data_to_nxdict, spice_scan_to_nxdict
from tavi.data.nexus_entry import NexusEntry
from tavi.data.spice_reader import _create_spicelogs


def test_nxdataset():
    name = NXdataset(ds="HFIR", type="NX_CHAR", EX_required="true")
    assert name.get_attr("type") == "NX_CHAR"
    assert name.get_attr("EX_required") == "true"
    assert name.get_dataset() == "HFIR"


def test_nxentry():
    nxsource = NXentry(
        name=NXdataset(ds="HFIR", type="NX_CHAR", EX_required="true"),
        probe=NXdataset(ds="neutron", type="NX_CHAR", EX_required="true"),
        NX_class="NXsource",
        EX_required="true",
    )
    assert len(nxsource.keys()) == 3
    assert type(nxsource["name"]) is NXdataset

    nxinstrument = NXentry(
        source=nxsource,
        NX_class="NXinstrument",
        EX_required="true",
    )
    assert nxinstrument["attrs"] == {"NX_class": "NXinstrument", "EX_required": "true"}
    assert nxinstrument["source"]["attrs"] == {"NX_class": "NXsource", "EX_required": "true"}
    assert nxinstrument["source"]["name"]["dataset"] == "HFIR"


def test_nxentry_from_daslogs():
    path_to_spice_data = "./test_data/exp424/Datafiles/CG4C_exp0424_scan0034.dat"
    spicelogs = _create_spicelogs(path_to_spice_data)
    nxdet = NXentry(
        data=NXdataset(ds=spicelogs["detector"], type="NX_INT", EX_required="true", unit="counts"),
        NX_class="NXdetector",
        EX_required="true",
    )
    assert nxdet["data"]["attrs"] == {"type": "NX_INT", "EX_required": "true", "unit": "counts"}


def test_add_dataset():
    path_to_spice_data = "./test_data/exp424/Datafiles/CG4C_exp0424_scan0034.dat"
    spicelogs = _create_spicelogs(path_to_spice_data)
    data = NXdataset(ds=spicelogs["detector"], type="NX_INT", EX_required="true", unit="counts")
    nxdet = NXentry(NX_class="NXdetector", EX_required="true")
    nxdet.add_dataset(key="data", ds=data)
    assert nxdet["data"]["attrs"] == {"type": "NX_INT", "EX_required": "true", "unit": "counts"}


def test_add_nonexisting_dataset():
    path_to_spice_data = "./test_data/exp424/Datafiles/CG4C_exp0424_scan0034.dat"
    spicelogs = _create_spicelogs(path_to_spice_data)

    nxdet = NXentry(
        none_data=NXdataset(ds=spicelogs.get("detector1"), type="NX_INT", EX_required="true", unit="counts"),
        NX_class="NXdetector",
        EX_required="true",
    )
    assert "none_data" not in nxdet.keys()

    nxdet = NXentry(NX_class="NXdetector", EX_required="true")
    none_data = NXdataset(ds=spicelogs.get("detector1"), type="NX_INT", EX_required="true", unit="counts")
    nxdet.add_dataset(key="none_data", ds=none_data)
    assert "none_data" not in nxdet.keys()


def test_spice_scan_to_nxdict():
    path_to_spice_data = "./test_data/exp424/Datafiles/CG4C_exp0424_scan0034.dat"
    nxdict = spice_scan_to_nxdict(path_to_spice_data)

    assert len(nxdict["SPICElogs"]) == 56
    assert nxdict["SPICElogs"]["attrs"]["scan"] == "34"
    assert nxdict["start_time"]["dataset"] == "2024-07-03T01:44:46-04:00"
    assert np.allclose(nxdict["instrument"]["monochromator"]["ei"]["dataset"][0:3], [4.9, 5, 5.1])

    entries = {"scan0034": nxdict}
    nexus_entries = {}
    for scan_num, scan_content in entries.items():
        content_list = []
        for key, val in scan_content.items():
            content_list.append((key, val))
        nexus_entries.update({scan_num: NexusEntry(content_list)})

    path_to_nexus_entry = "./test_data/spice_to_nxdict_test_scan0034.h5"
    for scan_num, nexus_entry in nexus_entries.items():
        nexus_entry.to_nexus(path_to_nexus_entry, scan_num)


def test_spice_data_to_nxdict():
    path_to_spice_data = "./test_data/exp424"
    nxdict = spice_data_to_nxdict(path_to_spice_data)

    assert len(nxdict) == 92
    assert nxdict["scan0034"]["SPICElogs"]["attrs"]["scan"] == "34"
