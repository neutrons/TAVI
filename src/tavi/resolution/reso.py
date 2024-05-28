import numpy as np
import numpy.linalg as la
from tavi.utilities import *
from tavi.tas import TAS


class Reso(TAS):
    """Resolution ellipoid calculator

    Attributs:
        STATUS (bool): True if resolution calculation is successful

    Methods:

    """

    def __init__(self):
        super().__init__()
        self.STATUS = False

    def ellipsoid_volume(self):
        """volume of the ellipsoid"""
        pass
