from typing import Literal, Optional

import numpy as np
from tavi.library.utilities import MotorAngles

class Goniometer:
    """
    Goniometer
    For Huber table, use type Y,-Z,X or Y,Z,-X
    For Four-Cricle in bisect mode, use type ?Y,Z,Y,bisect

    Attributes:
        type (str):
        sense (str): "+" if s2 is right-hand
        limits (dict): (min, max) of s2 angle

    Methods:
        _sense
        stacking_order
        mode (str): Mode for Four-Circle. "bisect" or "azimuthal"
        motor_senses: signs in the goniometer type string
        r_mat
        r_mat_inv
        set_limit
        validate_motor_positions

    """
    def __init__(self, type:str = "Y, -Z, X", sense: Optional[Literal["-", "+"]]= None):
        self.type = type
        self.sense = sense