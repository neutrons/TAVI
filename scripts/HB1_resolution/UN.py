import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axisartist import Axes

from tavi.data.scan import Scan
from tavi.instrument.tas import TAS
from tavi.plotter import Plot2D
from tavi.sample import Sample

ei = 13.500174
ef = 13.500173
instrument_config_json_path = "./test_data/IPTS31591_HB1_exp0917/hb1.json"
hb1 = TAS(fixed_ef=ef)
hb1.load_instrument_params_from_json(instrument_config_json_path)

path_to_spice_folder = "./test_data/exp974"
scan52 = Scan.from_spice(path_to_spice_folder, scan_num=52)
sample = Sample.from_scan(scan52)
hb1.mount_sample(sample)

#  -------------- check UB calculation ----------------
# [UBMode]
# Mode=1

# [Data]
# ScatteringPlaneVectors="1.000000,1.000000,0.000000,0.000000,0.000000,1.000000"
# Peak1="0.000000,0.000000,2.000000,-60.449408,-40.547250,3.000000,-1.500500,13.500174,13.500173"
# LatticeParams="4.890000,4.890000,4.890000,90.000000,90.000000,90.000000"
# Energy=13.500000

# [Matrices]
# UBMatrix="0.018308,0.033444,0.200913,-0.146030,-0.138468,0.036356,0.141985,-0.146725,0.011485"
# UBInverse="0.437778,-3.491889,3.395167,0.799708,-3.311050,-3.508492,4.804260,0.869342,0.274639"
# BMatrix="0.204499,-0.000000,-0.000000,0.000000,0.204499,-0.000000,0.000000,0.000000,0.204499"

# [AngleMode]
# Mode=0
# PlaneNormal="-0.052336,-0.026150,0.998287"
# InPlaneRef="-0.000000,-0.999657,-0.026186"
# UpperArc=-1.500500
# LowerArc=3.000000
# --------------------------------------------------------------
ub_mat_spice = sample.ub_conf.ub_mat

# angles1 = MotorAngles(two_theta=-60.449408, omega=-40.547250, sgl=3, sgu=-1.500500)
# peak1 = Peak(hkl=(0, 0, 2), angles=angles1)
# scattering_plane = ((1, 1, 0), (0, 0, 1))

# hb1.calculate_ub_matrix(peaks=peak1, scattering_plane=scattering_plane)
# print("UB from Scan:\n", ub_mat_spice)
# print("UB from calculation:\n", hb1.sample.ub_conf.ub_mat)

# ------------ genreate contour plot---------------
hkle_list = [(2, 2, ql, en) for ql in np.arange(0, 2, 0.2) for en in np.arange(0, 20, 3)]
# en_list = [e for e in np.arange(0, 50, 5)]
projection = ((1, 1, 0), (0, 0, 1), (1, -1, 0), "en")
rez_list = hb1.cooper_nathans(hkle=hkle_list, axes=projection)
p1 = Plot2D()
for rez in filter(None, rez_list):
    e_co = rez.get_ellipse(axes=(1, 3), PROJECTION=False)
    e_inco = rez.get_ellipse(axes=(1, 3), PROJECTION=True)
    p1.add_reso(e_co, c="k", linestyle="solid")
    p1.add_reso(e_inco, c="k", linestyle="dashed")

fig = plt.figure()
ax = fig.add_subplot(111, axes_class=Axes)

im = p1.plot(ax)
# ------------ genreate contour plot---------------
hkle_list = [(qh, qh, ql, 14) for ql in np.arange(-3, 3, 0.2) for qh in np.arange(-3, 3, 0.2)]
# en_list = [e for e in np.arange(0, 50, 5)]
projection = ((1, 1, 0), (0, 0, 1), (1, -1, 0), "en")
rez_list = hb1.cooper_nathans(hkle=hkle_list, axes=projection)
p2 = Plot2D()
for rez in filter(None, rez_list):
    e_co = rez.get_ellipse(axes=(0, 1), PROJECTION=False)
    e_inco = rez.get_ellipse(axes=(0, 1), PROJECTION=True)
    p2.add_reso(e_co, c="k", linestyle="solid")
    p2.add_reso(e_inco, c="k", linestyle="dashed")

fig = plt.figure()
ax = fig.add_subplot(111, axes_class=Axes)

im = p2.plot(ax)
# ---------- resolution for Q=(2,2,1), En=15 meV --------
rez0 = hb1.cooper_nathans(hkle=(2, 2, 1, 14), axes=projection)
# plotting
fig = plt.figure(figsize=(10, 6), constrained_layout=True)
rez0.plot_ellipses(fig)
plt.show()
