import numpy as np
from tavi.utilities import *
from tavi.sample.sample import Sample


class Xtal(Sample):
    """Singel crystal sample

    Attibutes:
        type (str): "xtal"

    """

    def __init__(self, lattice_params):
        super().__init__(lattice_params)
        self.type = "xtal"


if __name__ == "__main__":
    xtal = Xtal(lattice_params=(1, 1, 1, 90, 90, 120))

    print(xtal.b_mat() / 2 / np.pi)
    print(xtal.b_mat() @ np.array([0, 1, 0]) / 2 / np.pi)
