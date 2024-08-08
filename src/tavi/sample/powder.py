from tavi.sample.sample import Sample


# TODO
class Powder(Sample):
    """Powder sample

    Attibutes:
        type (str): "powder"

    """

    def __init__(self, lattice_params):
        super().__init__(lattice_params)
        self.type = "powder"

    @classmethod
    def from_json(cls, sample_params):
        """Alternate constructor from json"""
        lattice_params = (
            sample_params["a"],
            sample_params["b"],
            sample_params["c"],
            sample_params["alpha"],
            sample_params["beta"],
            sample_params["gamma"],
        )

        sample = cls(lattice_params=lattice_params)

        return sample
