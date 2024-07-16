import numpy as np
from tavi.sample.sample import Sample


class Powder(Sample):
    """Powder sample

    Attibutes:
        type (str): "powder"

    """

    def __init__(self, lattice_params):
        super().__init__(lattice_params)
        self.type = "powder"


if __name__ == "__main__":
    powder = Powder(lattice_params=(1, 1, 1, 90, 90, 120))

    print(powder.b_mat() / 2 / np.pi)
    print(powder.b_mat() @ np.array([0, 1, 0]) / 2 / np.pi)
