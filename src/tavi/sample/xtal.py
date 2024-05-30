import numpy as np
from tavi.utilities import *
from tavi.sample.sample import Sample


class Xtal(Sample):
    """Singel crystal sample

    Attibutes:
        type (str): "xtal"


    Methods:


    """

    def __init__(self, lattice_params):
        super().__init__(lattice_params)
        self.type = "xtal"
        self.ub_matrix = None

        self.i_star, self.j_star, self.k_star = self.reciprocal_basis()

    def set_UB(self, u=[1, 0, 0], v=[0, 0, 1]):
        self.u = u
        self.v = v


if __name__ == "__main__":
    xtal = Xtal(lattice_params=(3.574942, 3.574942, 5.663212, 90, 90, 120))

    print(xtal.b_mat() / 2 / np.pi)
    print(xtal.b_mat() @ np.array([0, 1, 0]) / 2 / np.pi)
