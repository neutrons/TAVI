from typing import Literal, Optional

import numpy as np

from tavi.instrument.components.component_base import TASComponent
from tavi.utilities import MotorAngles


class Goniometer(TASComponent):
    """
    Goniometer
    For Huber table, use type Y-ZX or YZ-X
    For Four-Cricle, use type ???
    """

    def __init__(
        self,
        param_dict: Optional[dict] = None,
        component_name: str = "goniometer",
    ):
        self.type: str = "Y-ZX"  # Y-mZ-X for (s1, sgl, sgu), Huber table at HB1A and HB3,
        self.sense: Literal["-", "+"] = "-"  # determines the sign of s2
        self.omega_limit: Optional[tuple]
        self.sgl_limit: Optional[tuple] = (-10, 10)
        self.sgu_limit: Optional[tuple] = (-10, 10)
        self.chi_limit: Optional[tuple]
        self.phi_limit: Optional[tuple]

        super().__init__(param_dict)
        self.component_name = component_name

    @property
    def _sense(self):
        match self.sense:
            case "+":
                return +1
            case "-":
                return -1
            case _:
                raise ValueError(f"sense {self.sense} needs to be either '+' or '-'.")

    @staticmethod
    def rot_x(nu):
        """rotation matrix about y-axis by angle nu

        Args:
            nu (float): angle in degrees

        Note:
            Using Mantid convention, beam along z, y is up, x in plane
        """

        angle = np.deg2rad(nu)
        c = np.cos(angle)
        s = np.sin(angle)
        mat = np.array(
            [
                [1, 0, 0],
                [0, c, -s],
                [0, s, c],
            ]
        )
        return mat

    @staticmethod
    def rot_y(omega):
        """rotation matrix about y-axis by angle omega

        Args:
            omega (float): angle in degrees

        Note:
            Using Mantid convention, beam along z, y is up, x in plane
        """

        angle = np.deg2rad(omega)
        c = np.cos(angle)
        s = np.sin(angle)
        mat = np.array(
            [
                [c, 0, s],
                [0, 1, 0],
                [-s, 0, c],
            ]
        )
        return mat

    @staticmethod
    def rot_z(mu):
        """rotation matrix about z-axis by angle mu

        Args:
            mu (float): angle in degrees

        Note:
            Using Mantid convention, beam along z, y is up, x in plane
        """

        angle = np.deg2rad(mu)
        c = np.cos(angle)
        s = np.sin(angle)
        mat = np.array(
            [
                [c, -s, 0],
                [s, c, 0],
                [0, 0, 1],
            ]
        )
        return mat

    def r_mat(self, angles: MotorAngles) -> np.ndarray:
        """Goniometer rotation matrix R

        Note:
            type refers to the axis of rotation for each motor,
            + for CCW, - for CW.
            The Cartesian coordinate Z is along the incoming beam,
            X in plane, Y up, right handed.
        """

        match self.type:  # (s1, sgl, sgu, chi, phi)
            case "Y-ZX":  # HB3
                r_mat = np.matmul(
                    Goniometer.rot_y(angles.omega),
                    np.matmul(
                        Goniometer.rot_z(-1 * angles.sgl),
                        Goniometer.rot_x(angles.sgu),
                    ),
                )
            case "YZ-X":  # CG4C ??
                r_mat = np.matmul(
                    Goniometer.rot_y(angles.omega),
                    np.matmul(
                        Goniometer.rot_z(angles.sgl),
                        Goniometer.rot_x(-1.0 * angles.sgu),
                    ),
                )
            case "YZY_bisect":
                r_mat = np.matmul(
                    Goniometer.rot_y(angles.omega),
                    np.matmul(
                        Goniometer.rot_z(angles.chi),
                        Goniometer.rot_y(angles.phi),
                    ),
                )
            case _:
                raise ValueError("Unknow goniometer type. Curruntly support Y-ZX and YZ-X")

        return r_mat

    def r_mat_inv(
        self,
        angles: MotorAngles,
    ) -> np.ndarray:
        """inverse of rotation matrix, equivalent to transpose since R is unitary"""
        # return np.linalg.inv(self.r_mat(angles))
        return self.r_mat(angles).T

    def angles_from_r_mat(
        self,
        r_mat: np.ndarray,
    ) -> MotorAngles:
        """Calculate goniometer angles from the R matrix

        Note:
            range of np.arcsin is -pi/2 to pi/2
            range of np.atan2 is -pi to pi
        """

        def stacking_order_yzx(r_mat):
            """Huber table, return angles in degrees"""
            # sgl1 = np.arcsin(r_mat[1, 0]) * rad2deg
            # sgl2 = np.arccos(np.sqrt(r_mat[0, 0] ** 2 + r_mat[2, 0] ** 2)) * rad2deg
            sgl_rad = np.arctan2(r_mat[1, 0], np.sqrt(r_mat[0, 0] ** 2 + r_mat[2, 0] ** 2))
            sgl = np.rad2deg(sgl_rad)

            # sgu1 = np.arcsin(-r_mat[1, 2] / np.sqrt(r_mat[0, 0] ** 2 + r_mat[2, 0] ** 2)) * rad2deg
            # sgu2 = np.arccos(r_mat[1, 1] / np.sqrt(r_mat[0, 0] ** 2 + r_mat[2, 0] ** 2)) * rad2deg
            sgu_rad = np.arctan2(
                -r_mat[1, 2] / np.sqrt(r_mat[0, 0] ** 2 + r_mat[2, 0] ** 2),
                r_mat[1, 1] / np.sqrt(r_mat[0, 0] ** 2 + r_mat[2, 0] ** 2),
            )
            sgu = np.rad2deg(sgu_rad)

            # omega1 = np.arcsin(-r_mat[2, 0] / np.sqrt(r_mat[0, 0] ** 2 + r_mat[2, 0] ** 2)) * rad2deg
            # omega2 = np.arccos(r_mat[0, 0] / np.sqrt(r_mat[0, 0] ** 2 + r_mat[2, 0] ** 2)) * rad2deg
            omega_rad = np.arctan2(
                -r_mat[2, 0] / np.sqrt(r_mat[0, 0] ** 2 + r_mat[2, 0] ** 2),
                r_mat[0, 0] / np.sqrt(r_mat[0, 0] ** 2 + r_mat[2, 0] ** 2),
            )
            omega = np.rad2deg(omega_rad)
            return omega, sgl, sgu

        # TODO four-circle is chi circle parallel to beam
        def stacking_order_yzy(r_mat):

            pass

        # TODO four-circle is chi circle perpendicular to beam
        def stacking_order_yxy(r_mat):
            pass

        match self.type:
            case "Y-ZX":  # Y-mZ-X (s1, sgl, sgu) for HB1A and HB3,
                omega, sgl, sgu = stacking_order_yzx(r_mat)
                angles = MotorAngles(None, omega, -1 * sgl, sgu, None, None)
                self.validate_motor_positions(angles)

            case "YZ-X":  # Y-Z-mX (s1, sgl, sgu) for CG4C
                omega, sgl, sgu = stacking_order_yzx(r_mat)
                angles = MotorAngles(None, omega, sgl, -1 * sgu, None, None)
                self.validate_motor_positions(angles)

            case "YZ-Y_bisect":
                omega, chi, phi = stacking_order_yzy(r_mat)
                angles = MotorAngles(None, omega, None, None, chi, phi)
                self.validate_motor_positions(angles)

            case _:
                raise ValueError("Unknow goniometer type.  Curruntly support Y-ZX and YZ-X.")

        return angles

    def set_limit(
        self,
        motor_name: Literal["omega", "sgl", "sgu", "chi", "phi"],
        motor_range: tuple[float, float] = (-180, 180),
    ):
        "set goiometer motor limt"
        setattr(self, motor_name + "_limit", motor_range)

    # TODO
    def validate_motor_positions(self, angles: MotorAngles):
        "check if all goiometer motors are within the limits"
        pass
