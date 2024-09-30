from tavi.data.nxdict import NXdataset, NXentry
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
    spicelogs = _create_spicelogs(path_to_spice_data)
