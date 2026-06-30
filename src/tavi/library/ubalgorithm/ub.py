"""Calculating ub matrix related algorithms."""

import numpy as np

from tavi.library.experiment.peak import DataPoint
from tavi.library.experiment.utilities import q_lab
from tavi.library.geometry.sample import Sample
from tavi.library.Instrument.instrument import Instrument


class UBAlgorithm:
    """Calculate UB."""

    def __init__(self, sample: Sample, instrument: Instrument) -> None:
        """Init ub algorithm."""
        self.sample = sample
        self.instrument = instrument

        # -----------Calculate UB Matrix-----------

    def find_u_from_one_peak_and_scattering_plane(
        self,
        peak: DataPoint,
        scattering_plane: tuple[tuple, tuple],
        ei: float,
        ef: float,
    ) -> np.ndarray:
        """Calculate U matrix from one peak and a scattering plane."""
        B = self.sample.ol.B
        t1_c = B @ peak.hkl
        vector1, vector2 = scattering_plane

        coeff1, coeff2 = vector1 @ peak.hkl, vector2 @ peak.hkl
        # construct t2_c as orthogonal as possible from t1_c
        if np.abs(coeff1) > np.abs(coeff2):  # vector1 is closer to t1_c, use vector2
            t3_c = np.cross(B @ peak.hkl, B @ vector2 * np.sign(coeff1))  # the sign guarantee right-hand-rule
        else:
            t3_c = np.cross(B @ peak.hkl, B @ vector1 * np.sign(coeff2))
        t2_c = np.cross(t3_c, t1_c)

        T_c = np.array(
            [
                t1_c / np.linalg.norm(t1_c),
                t2_c / np.linalg.norm(t2_c),
                t3_c / np.linalg.norm(t3_c),
            ]
        ).T

        q_lab1 = q_lab(ei, ef, peak.angles.angles_dict["two_theta"])
        t1_v = np.linalg.inv(self.instrument.goni.r_mat(peak.angles)) @ q_lab1
        # Rotate to a mantid coordinate system.
        t3_v = np.linalg.inv(self.instrument.goni.r_mat(peak.angles)) @ np.array([0, 1, 0]).T
        t2_v = np.cross(t3_v, t1_v)
        T_v = np.array(
            [
                t1_v / np.linalg.norm(t1_v),
                t2_v / np.linalg.norm(t2_v),
                t3_v / np.linalg.norm(t3_v),
            ]
        ).T

        u_mat = T_v @ T_c.T
        return u_mat

    def find_u_from_two_peaks(self, peaks: tuple[DataPoint, DataPoint]) -> np.ndarray:
        """
        Calculate U matrix from two peaks.

        r_mat can be removed later when goniometer is implemented as it can be calculated from
        peak.angles.

        Follow Eq.76-81 and Eq.83-88. We assume q_3 is perpendicular from the two peaks
        """
        peak1, peak2 = peaks
        B = self.sample.ol.B

        t1_c = B @ peak1.hkl
        t3_c = np.cross(B @ peak1.hkl, B @ peak2.hkl)
        t2_c = np.cross(t3_c, t1_c)

        T_c = np.array(
            [
                t1_c / np.linalg.norm(t1_c),
                t2_c / np.linalg.norm(t2_c),
                t3_c / np.linalg.norm(t3_c),
            ]
        ).T

        # We need to create vectors t1_v, t2_v, t3_v
        q_lab_1 = q_lab(peak1.ei, peak1.ef, peak1.angles.angles_dict["two_theta"])
        q_lab_2 = q_lab(peak1.ei, peak2.ef, peak2.angles.angles_dict["two_theta"])

        # In identical fashion as described above Eq.79
        t1_v = np.linalg.inv(self.instrument.goni.r_mat(peak1.angles)) @ q_lab_1 / (2 * np.pi)
        t3_v = np.cross(
            np.linalg.inv(self.instrument.goni.r_mat(peak1.angles)) @ q_lab_1 / (2 * np.pi),
            np.linalg.inv(self.instrument.goni.r_mat(peak2.angles)) @ q_lab_2 / (2 * np.pi),
        )
        t2_v = np.cross(t3_v, t1_v)
        T_v = np.array(
            [
                t1_v / np.linalg.norm(t1_v),
                t2_v / np.linalg.norm(t2_v),
                t3_v / np.linalg.norm(t3_v),
            ]
        ).T

        u_mat = T_v @ T_c.T
        return u_mat

    def find_ub_from_three_peaks(self, peaks: tuple[DataPoint, DataPoint, DataPoint]) -> np.ndarray:
        """
        Calculate U matrix from three peaks.

        r_mat can be removed later when goniometer is implemented as it can be calculated from
        peak.angles.

        Follow Eq.83-90.
        """
        peak1, peak2, peak3 = peaks

        r_mat = self.instrument.goni.r_mat
        # we directly use the three peaks as t1_c, t2_c and t3_c
        V = np.array([peak1.hkl, peak2.hkl, peak3.hkl]).T
        q_lab_1 = q_lab(peak1.ei, peak1.ef, peak1.angles.angles_dict["two_theta"])
        q_lab_2 = q_lab(peak2.ei, peak2.ef, peak2.angles.angles_dict["two_theta"])
        q_lab_3 = q_lab(peak3.ei, peak3.ef, peak3.angles.angles_dict["two_theta"])

        q1_v = np.linalg.inv(r_mat(peak1.angles)) @ q_lab_1 / (2 * np.pi)
        q2_v = np.linalg.inv(r_mat(peak2.angles)) @ q_lab_2 / (2 * np.pi)
        q3_v = np.linalg.inv(r_mat(peak3.angles)) @ q_lab_3 / (2 * np.pi)
        Q_v = np.array([q1_v, q2_v, q3_v]).T

        ub_mat = Q_v @ np.linalg.inv(V)
        return ub_mat

    def find_ub_from_multiple_peaks(self, peaks: tuple[DataPoint, ...]) -> np.ndarray:
        """
        Calculate U matrix from three peaks.

        r_mat can be removed later when goniometer is implemented as it can be calculated from
        peak.angles.

        Follow Eq.89-98.
        """
        n = len(peaks)
        Q_v = np.zeros((3, 3))
        VV = np.zeros((3, 3))
        r_mat = self.instrument.goni.r_mat

        for i in range(n):
            hkl = peaks[i].hkl
            q_lab_i = q_lab(peaks[i].ei, peaks[i].ef, peaks[i].angles.angles_dict["two_theta"])
            q_v_i = np.linalg.inv(r_mat(peaks[i].angles)) @ q_lab_i / (2 * np.pi)
            for j in range(3):
                for k in range(3):
                    Q_v[j, k] += q_v_i[k] * hkl[j]
                    VV[j, k] += hkl[k] * hkl[j]
        ub_mat = Q_v.T @ np.linalg.inv(VV).T
        return ub_mat

    def plane_normal_from_two_peaks(self, peaks: tuple[DataPoint, DataPoint]) -> tuple[np.ndarray, np.ndarray]:
        """
        Calculate plane_normal and in_plane reflection.

        Both are vectors representing peaks in Q_lab.
        """
        peak1, peak2 = peaks
        u_mat = self.sample.ol.U
        b_mat = self.sample.ol.B
        t1_c = b_mat @ peak1.hkl
        t3_c = np.cross(b_mat @ peak1.hkl, b_mat @ peak2.hkl)

        # Eq. 79(3) calculate t3_v from U and t3_c
        plane_normal = u_mat @ t3_c / np.linalg.norm(t3_c)
        # if y is pointing down, set it to point up
        plane_normal = -plane_normal if plane_normal[1] < 0 else plane_normal
        in_plane_ref = u_mat @ t1_c / np.linalg.norm(t1_c)
        return (plane_normal, in_plane_ref)
