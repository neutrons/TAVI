import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axisartist import Axes

from tavi.data.nexus_builder import NXdataset, NXentry
from tavi.data.nexus_entry import NexusEntry
from tavi.data.tavi import TAVI
from tavi.instrument.resolution.cooper_nathans_bak import CooperNathans
from tavi.plotter import Plot2D
from tavi.sample.xtal import Xtal

instrument_config_json_path = "./src/tavi/instrument/instrument_params/cg4c.json"
tas = CooperNathans(SPICE_CONVENTION=False)
tas.load_instrument_params_from_json(instrument_config_json_path)

sample_json_path = "./test_data/test_samples/nitio3.json"
sample = Xtal.from_json(sample_json_path)
tas.mount_sample(sample)


tavi = TAVI()

# calculate resolution ellipses
R0 = False
hkl_list = [(qh, 0, 0) for qh in np.arange(1, 2.01, 0.1)]
ef = 4.8
ei_list = [e + ef for e in np.arange(0, 5.1, 0.2)]

# genreate plot
p = Plot2D()

for i, hkl in enumerate(hkl_list, 1):
    rez_list = tas.rez(hkl_list=hkl, ei=ei_list, ef=ef, R0=R0)
    sz = len(rez_list)
    rez_entry = NXentry(
        hkl=NXdataset(ds=[rez.hkl for rez in rez_list]),
        en=NXdataset(ds=[rez.en for rez in rez_list]),
        rez_mat=NXdataset(ds=[rez.mat for rez in rez_list]),
        rez_r0=NXdataset(ds=[rez.r0 for rez in rez_list]),
    )
    tavi.processed_data.update(NexusEntry._dict_to_nexus_entry({f"scan{i:04}": rez_entry}))

    for rez in rez_list:
        e_co = rez.get_ellipse(axes=(0, 3), PROJECTION=False)
        e_inco = rez.get_ellipse(axes=(0, 3), PROJECTION=True)
        if e_co is not None:
            p.add_reso(e_co, c="k", linestyle="solid")
        if e_inco is not None:
            p.add_reso(e_inco, c="k", linestyle="dashed")

tavi.save("./test_data/rez_dispH.h5")


fig = plt.figure()
ax = fig.add_subplot(111, axes_class=Axes)
im = p.plot(ax)
plt.show()
