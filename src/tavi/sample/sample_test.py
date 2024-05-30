from tavi.utilities import *
from tavi.sample.xtal import Xtal

test_xtal = Xtal(lattice_params=(5.3995, 5.64, 11.75, 90, 90, 90))

test_xtal.shape = "cylindrical"
test_xtal.width = 1.0 * cm2A
test_xtal.height = 1.0 * cm2A
test_xtal.depth = 1.0 * cm2A
test_xtal.mosaic = 30 * min2rad  # horizontal mosaic
test_xtal.mosaic_v = 30 * min2rad  # vertical mosaic
test_xtal.ub_matrix = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
