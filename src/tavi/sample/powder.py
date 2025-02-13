from tavi.sample import Sample


# TODO
class Powder(Sample):
    """Powder sample

    Attibutes:
        type (str): "powder"

    """

    def __init__(self, lattice_params):
        super().__init__(lattice_params)
        self.type = "powder"
