import numpy as np

from tavi.sample.xtal import Xtal

# from tavi.utilities import *

# test_xtal = Xtal(lattice_params=(5.3995, 5.64, 11.75, 90, 90, 90))
test_xtal = Xtal(lattice_params=(3.574924, 3.574924, 5.663212, 90, 90, 120))
# test_xtal = Xtal(lattice_params=(5.663212, 5.663212, 5.663212, 90, 90, 90))

test_xtal.shape = "cylindrical"
test_xtal.width = 1.0  # * cm2angstrom
test_xtal.height = 1.0  # * cm2angstrom
test_xtal.depth = 1.0  # * cm2angstrom
test_xtal.mosaic_h = 30  # * min2rad  # horizontal mosaic
test_xtal.mosaic_v = 30  # * min2rad  # vertical mosaic

ub_matrix_spice = np.array(
    [
        [0.053821, 0.107638, 0.166485],
        [-0.164330, -0.304247, 0.058788],
        [0.272815, -0.013290, 0.002566],
    ]
)
# ub matrix in Mantid convention
test_xtal.ub_matrix = np.array(
    [
        ub_matrix_spice[0],
        ub_matrix_spice[2],
        -ub_matrix_spice[1],
    ]
)

test_xtal.plane_normal_spice = np.array([0.000009, -0.043637, 0.999047])
test_xtal.plane_normal = np.array(
    [
        test_xtal.plane_normal_spice[0],
        test_xtal.plane_normal_spice[2],
        -test_xtal.plane_normal_spice[1],
    ]
)
