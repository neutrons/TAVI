import xml.etree.ElementTree as ET

import numpy as np

import re

_SHAPE = re.compile(r"^(\w+?)(?:\[([\d,]+)\])?$")


_CASTS = {"FLOAT32": float, "FLOAT64": float, "INT32": int, "INT64": int}


def _value(elem: ET.Element):
    """Convert one SPICE XML leaf to a scalar, array, or image per its `type` attribute."""
    text = (elem.text or "").strip()
    match = _SHAPE.match(elem.get("type", ""))
    cast = _CASTS.get(match.group(1)) if match else None
    if cast is None:  # untyped: Header strings, DetectorDataFilePath
        return text
    if match.group(2) is None:
        return cast(text)
    shape = tuple(int(n) for n in match.group(2).split(","))
    return np.fromstring(text, dtype=np.int64 if cast is int else np.float64, sep="\t").reshape(shape)


def _walk(elem: ET.Element):
    """Recurse into section elements; leaves become typed values."""
    if len(elem) == 0:
        return _value(elem)
    return {child.tag: _walk(child) for child in elem}


def parse_spice_xml(path: str) -> dict:
    """Parse one DEMAND/HB3A SPICE XML point file into nested dicts of typed values."""
    root = ET.parse(path).getroot()
    return {"file": dict(root.attrib)} | {child.tag: _walk(child) for child in root}

data = parse_spice_xml('/Users/qmc/Neutrons/tasdata/DEMAND/IPTS-32174/exp982/Datafiles/HB3A_exp982_scan0063_0191.xml')
print(
    len(data['Data']['Detector'][0])
    )
camera_matrix = data['Data']['Detector'].T
from matplotlib import pyplot as plt
plt.imshow(camera_matrix, cmap = 'viridis', origin = 'lower')
plt.colorbar()
# plt.show()

wavelength = 1.533
k = 2*np.pi/wavelength

theta = -30 # in degrees
theta = np.radians(theta)
q_l = [-k * np.sin(theta), 0, k - k * np.cos(theta)]
ub = np.array([
    [0.170739, 0.085370, 0.030710], 
    [-0.203490, -0.101745, 0.025767],
    [0.000000, -0.230044, 0.000000]
])
from tavi.library.experiment.utilities import spice_to_mantid, rot_z,rot_y
ub = spice_to_mantid(ub, version="old")
phi = -39.9935
chi=-0.0031
omega = -147.4995

r_mat = (
    rot_y(omega)
    @ rot_z(chi)
    @ rot_y(phi)
)
hkl = np.linalg.inv(ub) @ np.linalg.inv(r_mat) @ q_l / (2*np.pi)
print(hkl)
#-0.9373    -0.0001    -5.6899

def rot_theta_y(xyz: tuple, theta):
    """theta is mantid FIG.1, sometimes indicated as 2theta, in degrees."""
    theta *= -1
    theta = np.radians(theta)
    xyz = np.asarray(xyz)
    new_xyz = rot_y(theta) @ xyz.T
    return new_xyz.T

def translate(xyz, theta, d):
    """theta is mantid FIG.1, sometimes indicated as 2theta, in degrees."""
    theta *= -1
    theta = np.radians(theta)
    x,y,z = xyz
    z += d * np.cos(np.abs(theta))
    x += d *np.sin(theta)
    return np.array([x,y,z])

def get_polar_coordinate(xyz):
    """Get theta and phi"""
    x,y,z = xyz
    return np.degrees(np.arctan(x/z))

xyz = [0, 0, 0]
new_xyz = translate(rot_theta_y(xyz, 30), 30, 409.9628)

print(get_polar_coordinate(new_xyz))
