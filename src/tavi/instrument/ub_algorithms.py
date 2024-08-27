# -*- coding: utf-8 -*-

import numpy as np


def find_ub(self, peaks, angles, ei=13.5, ef=None):
    """calculate UB matrix from peaks and motor positions

    Args:
        peaks (list)
        angles (list)
        ei (float): incident neutron energy, in meV
        ef (float): final neutron energy, in meV

    """

    self.sample.ub_peaks = peaks
    self.sample.ub_angles = angles

    ki = np.sqrt(ei / ksq2eng)
    if ef is None:
        kf = ki
    else:
        kf = np.sqrt(ef / ksq2eng)

    if not len(peaks) == len(angles):
        print("Number of peaks and angles provided do not match.")

    if len(peaks) == 2:
        b_mat = self.sample.b_mat()
        (
            u_mat,
            in_plane_ref,
            plane_normal,
        ) = self._find_u_from_2peaks(peaks, angles, ki, kf)
        ub_matrix = u_mat @ b_mat

        self.sample.in_plane_ref = in_plane_ref
        self.sample.plane_normal = plane_normal

    # TODO
    elif len(peaks) == 3:  # find_ub_from_3peaks
        pass
    # TODO
    elif len(peaks) > 3:  # find_ub_from_mulitple_peaks
        pass
    else:
        print("I don't even know what you're doing.")

    self.sample.ub_matrix = ub_matrix
    inv_ub_matrix = np.linalg.inv(ub_matrix)
    self.sample.inv_ub_matrix = inv_ub_matrix

    # print(np.round(ub_matrix, 6))


def find_angles(self, peak, ei=13.5, ef=None):
    """calculate motor positions for a given peak if UB matrix has been determined

    Args:
        peak (tuple): (h, k, l) of a peak
        ei (float): incident neutron energy, in meV
        ef (float): final neutron energy, in meV

    """
    S2_MIN_DEG = 1

    hkl = np.array(peak)
    ki = np.sqrt(ei / ksq2eng)
    if ef is None:
        kf = ki
    else:
        kf = np.sqrt(ef / ksq2eng)

    b_mat = self.sample.b_mat()

    q_sq = 4 * np.pi**2 * hkl.T @ b_mat.T @ b_mat @ hkl
    q_norm = np.sqrt(q_sq)

    # two_theta = np.arccos((ki**2 + kf**2 - q_sq) / (2 * ki * kf)) * self.goniometer.sense
    two_theta = get_angle(ki, kf, q_norm)

    if two_theta is None:
        print(f"Triangle cannot be closed at q={hkl}, en={ei-ef} meV.")
        angles = None
    elif np.rad2deg(two_theta) < S2_MIN_DEG:
        print(f"s2 is smaller than {S2_MIN_DEG} deg at q={hkl}.")
        angles = None
    else:
        two_theta = np.rad2deg(two_theta) * self.goniometer.sense
        q = self.sample.ub_matrix @ hkl
        t1 = q / np.linalg.norm(q)

        eps = 1e-8  # zero
        plane_normal = self.sample.plane_normal
        in_plane_ref = self.sample.in_plane_ref

        # Minimal tilts
        if np.dot(t1, plane_normal) < eps:  # t1 in plane
            t3 = plane_normal
            t2 = np.cross(t3, t1)
        else:  # t1 not in plane, need to change tilts
            if np.linalg.norm(np.cross(plane_normal, t1)) < eps:
                # oops, t1 along plane_normal
                t2 = in_plane_ref
                t3 = np.cross(t1, t2)
            else:
                t2p = np.cross(plane_normal, t1)
                t3 = np.cross(t1, t2p)
                t2 = np.cross(t3, t1)

        t_mat = np.array(
            [t1, t2 / np.linalg.norm(t2), t3 / np.linalg.norm(t3)],
        ).T

        t_mat_inv = np.linalg.inv(t_mat)

        q_lab1 = TripleAxisSpectrometer.q_lab(two_theta, ki, kf) / q_norm
        q_lab2 = np.array([q_lab1[2], 0, -q_lab1[0]])
        q_lab3 = np.array([0, 1, 0])

        q_lab_mat = np.array([q_lab1, q_lab2, q_lab3]).T
        r_mat = q_lab_mat @ t_mat_inv

        angles = np.round((two_theta,) + self.goniometer.angles_from_r_mat(r_mat), 6)

    return angles
