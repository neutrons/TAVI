import numpy as np

from tavi.sample.xtal import Xtal
from tavi.utilities import cm2angstrom, min2rad

# test_xtal = Xtal(lattice_params=(5.3995, 5.64, 11.75, 90, 90, 90))
nitio3 = Xtal(lattice_params=(5.034785, 5.034785, 13.812004, 90, 90, 120))
# test_xtal = Xtal(lattice_params=(5.663212, 5.663212, 5.663212, 90, 90, 90))

nitio3.shape = "cylindrical"
nitio3.width = 1.0 * cm2angstrom
nitio3.height = 1.0 * cm2angstrom
nitio3.depth = 1.0 * cm2angstrom
nitio3.mosaic_h = 30 * min2rad  # horizontal mosaic
nitio3.mosaic_v = 30 * min2rad  # vertical mosaic

ub_matrix_spice = np.array(
    [
        [-0.016965, -0.026212, -0.071913],
        [-0.201388, -0.193307, 0.007769],
        [-0.108415, 0.120600, -0.003178],
    ]
)
# ub matrix in Mantid convention
nitio3.ub_matrix = np.array(
    [
        ub_matrix_spice[0],
        ub_matrix_spice[2],
        -ub_matrix_spice[1],
    ]
)
nitio3.plane_normal_spice = np.array([-0.04032, 0.035237, 0.998565])
nitio3.plane_normal = np.array(
    [
        nitio3.plane_normal_spice[0],
        nitio3.plane_normal_spice[2],
        -nitio3.plane_normal_spice[1],
    ]
)
