from typing import Literal, Optional

import numpy as np

from tavi.instrument.tas_cmponents import TASComponent


class Goniometer(TASComponent):
    """Goniometer table, type = Y-ZX or YZ-X"""

    def __init__(
        self,
        param_dict: Optional[dict] = None,
        component_name: str = "goniometer",
    ):
        self.type: str = "Y-ZX"  # Y-mZ-X for Huber stage at HB1A and HB3
        self.sense: Literal[-1, +1] = -1

        super().__init__(param_dict)
        self.component_name = component_name

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

    def r_mat(
        self,
        angles_deg: tuple[float, float, float],
    ):
        "Goniometer rotation matrix R"

        omega, sgl, sgu = angles_deg  # s2, s1, sgl, sgu
        match self.type:
            case "Y-ZX":  # HB3
                r_mat = Goniometer.rot_y(omega) @ Goniometer.rot_z(-1 * sgl) @ Goniometer.rot_x(sgu)
            case "YZ-X":  # CG4C ??
                r_mat = Goniometer.rot_y(omega) @ Goniometer.rot_z(sgl) @ Goniometer.rot_x(-1 * sgu)
            case _:
                r_mat = None
                print("Unknow goniometer type. Curruntly support Y-ZX and YZ-X")

        return r_mat

    def r_mat_inv(
        self,
        angles: tuple[float, float, float],
    ):
        """inverse of rotation matrix"""
        # return np.linalg.inv(self.r_mat(angles))
        return self.r_mat(angles).T

    def angles_from_r_mat(self, r_mat):
        """Calculate goniometer angles from the R matrix

        Note:
            range of np.arcsin is -pi/2 to pi/2
            range of np.atan2 is -pi to pi
        """

        match self.type:
            case "Y-ZX" | "YZ-X":  # Y-mZ-X (s1, sgl, sgu) for HB1A and HB3, Y-Z-mX (s1, sgl, sgu) for CG4C
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

                match self.type:
                    case "Y-ZX":
                        angles = (omega, -1 * sgl, sgu)
                    case "YZ-X":
                        angles = (omega, sgl, -1 * sgu)

            case _:
                angles = None
                print("Unknow goniometer type.  Curruntly support Y-ZX and YZ-X.")

        return angles
