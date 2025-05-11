import concurrent.futures
from time import time

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse


def quadric_proj(quadric, idx):
    """projects along one axis of the quadric"""

    # delete if orthogonal
    zero = 1e-8
    if np.abs(quadric[idx, idx]) < zero:
        return np.delete(np.delete(quadric, idx, axis=0), idx, axis=1)

    # row/column along which to perform the orthogonal projection
    vec = 0.5 * (quadric[idx, :] + quadric[:, idx])  # symmetrise if not symmetric
    vec /= np.sqrt(quadric[idx, idx])  # normalise to indexed component
    proj_op = np.outer(vec, vec)  # projection operator
    ortho_proj = quadric - proj_op  # projected quadric

    return np.delete(np.delete(ortho_proj, idx, axis=0), idx, axis=1)


def incoh_sigma(mat, axis):
    """Incoherent sigma"""
    idx = int(axis)

    for i in (3, 2, 1, 0):
        if not i == idx:
            mat_proj = quadric_proj(mat, i)

    return 1 / np.sqrt(np.abs(mat_proj[0, 0]))


def model_disp(vq1, vq2, vq3):
    """return energy for given Q points
    3d FM J=-1 meV S=1, en=6*S*J*(1-cos(Q))
    """
    sj = 1
    gamma_q = (np.cos(2 * np.pi * vq1) + np.cos(2 * np.pi * vq2) + np.cos(2 * np.pi * vq3)) / 3
    return 6 * sj * (1 - gamma_q)


def model_inten(vq1, vq2, vq3):
    """return intensity for given Q points
    3d FM J=-1 meV S=1, inten = S/2 for all Qs
    """
    return np.ones_like(vq1, dtype=float) / 2


def resolution_matrix(qx0, qy0, qz0, en0):
    """Fake resoltuion matrix mat and prefactor r0
    r0 is a constant, rez_mat is a symmatric positive 4 by 4 matrix
    """

    sz = np.shape(qx0)

    def rotation_matrix_4d(theta_deg):
        theta = np.radians(theta_deg)
        c = np.cos(theta)
        s = np.sin(theta)
        return np.array(
            [
                [c, 0, 0, -s],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
                [s, 0, 0, c],
            ]
        )

    sigma1, sigma2 = 0.3, 0.01
    sigma3 = sigma4 = 1e-6
    angle = -80
    mat = np.array(
        [
            [1 / sigma1**2, 0, 0, 0],
            [0, 1 / sigma3**2, 0, 0],
            [0, 0, 1 / sigma4**2, 0],
            [0, 0, 0, 1 / sigma2**2],
        ]
    )
    rez_mat = rotation_matrix_4d(angle).T @ mat @ rotation_matrix_4d(angle)
    r0 = 1
    return np.broadcast_to(r0, sz), np.broadcast_to(rez_mat, sz + (4, 4))


def plot_rez_ellipses(ax):
    sigma1, sigma2 = 0.3, 0.02
    angle = 80
    for i in range(3):

        ax.add_artist(
            Ellipse(
                xy=(0, 0),
                width=sigma1 * 2 * (i + 1),
                height=sigma2 * 2 * (i + 1),
                angle=angle,
                edgecolor="w",
                facecolor="none",
                label=f"{i+1}-sigma",
            )
        )


# @njit(parallel=True)
def convolutioin(qh, qk, ql, en):

    # n_measrement = np.shape(qe_mesh)[1]
    # inten_measure = []

    # for i in range(n_measrement):
    # qh, qk, ql, en = qe_mesh[:, :]
    # ----------------------------------------------------
    # calculate resolution matrix for all points
    # ----------------------------------------------------
    r0, mat = resolution_matrix(qh, qk, ql, en)
    # ----------------------------------------------------
    # calculate the incoherent sigmas for three Q directions
    # ----------------------------------------------------
    sigma_qh = incoh_sigma(mat, 0)
    sigma_qk = incoh_sigma(mat, 1)
    sigma_ql = incoh_sigma(mat, 2)

    min_qh, max_qh = qh - 1.5 * sigma_qh, qh + 1.5 * sigma_qh
    min_qk, max_qk = qk - 1.5 * sigma_qk, qk + 1.5 * sigma_qk
    min_ql, max_ql = ql - 1.5 * sigma_ql, ql + 1.5 * sigma_ql

    # create 3D mesh
    pts_q = 20
    step_qh = (max_qh - min_qh) / pts_q
    step_qk = (max_qk - min_qk) / pts_q
    step_ql = (max_ql - min_ql) / pts_q

    list_qh = np.linspace(min_qh, max_qh, pts_q)
    list_qk = np.linspace(min_qk, max_qk, pts_q)
    list_ql = np.linspace(min_ql, max_ql, pts_q)

    mesh_q = np.meshgrid(list_qh, list_qk, list_ql, indexing="ij")
    vqh, vqk, vql = np.array([np.ravel(v) for v in mesh_q])
    # ----------------------------------------------------
    # calculate weight based on the dispersion
    # ----------------------------------------------------
    disp = model_disp(vqh, vqk, vql)  # size of pts_q^3
    # calculate weight from resolution function
    vqe = np.stack((vqh - qh, vqk - qk, vql - ql, disp - en), axis=1)  # size of (pts_q^3, 4)
    prod = np.einsum("ij,jk,ik->i", vqe, mat, vqe)
    weights = np.exp(-prod / 2)
    # ----------------------------------------------------
    # trim the corners
    # ----------------------------------------------------
    cut_off = 0.03
    idx = weights > np.max(weights) * cut_off

    if not np.any(idx):  # zero intensity
        # inten_measure.append(0.0)
        return 0.0

    else:
        vqe_filtered = vqe[idx]
        inten = model_inten(vqe_filtered[:, 0], vqe_filtered[:, 1], vqe_filtered[:, 2])
        weights_filtered = weights[idx]
        # normalization
        det = np.linalg.det(mat)
        inten_sum = np.sum(inten * weights_filtered) * step_qh * step_qk * step_ql
        # inten_measure.append(inten_sum * np.sqrt(det) / (2 * np.pi) ** 2 * r0)
        return inten_sum * np.sqrt(det) / (2 * np.pi) ** 2 * r0
    # return np.asarray(inten_measure)


if __name__ == "__main__":
    # ----------------------------------------------------
    # points being measured
    # qe_mesh has the dimension (4, n_pts_of_measurement)
    # flatten for meshed measurement
    # ----------------------------------------------------
    q1_min, q1_max, q1_step = -2, 2, 0.02
    en_min, en_max, en_step = -2, 6, 0.1
    q2 = 0
    q3 = 0

    q1 = np.linspace(q1_min, q1_max, int((q1_max - q1_min) / q1_step) + 1)
    en = np.linspace(en_min, en_max, int((en_max - en_min) / en_step) + 1)
    vq1, vq2, vq3, ven = np.meshgrid(q1, q2, q3, en, indexing="ij")

    sz = np.shape(vq1)  # sz = (n_q1, n_q2, n_q3 , n_en)
    qe_mesh = np.array([np.ravel(v) for v in (vq1, vq2, vq3, ven)])

    n_measrement = np.shape(qe_mesh)[1]

    t0 = time()
    num_worker = 4
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_worker) as executor:
        results = executor.map(convolutioin, *qe_mesh)
    print(f"Convolution completed in {(t1:=time())-t0:.4f} s")

    measurement_inten = np.array(list(results)).reshape(sz)
    # total intensity should be close to S/2 *(q1_max - q1_min) * 2p*i
    total_intent = np.sum(measurement_inten) * q1_step * en_step * 2 * np.pi

    # ----------------------------------------------------
    # plot 2D contour
    # ----------------------------------------------------
    fig, ax = plt.subplots()
    idx = np.s_[:, 0, 0, :]
    img = ax.pcolormesh(vq1[idx], ven[idx], measurement_inten[idx], cmap="turbo", vmin=0, vmax=0.5)

    ax.grid(alpha=0.6)
    ax.set_xlabel("Q1")
    ax.set_ylabel("En")
    ax.set_xlim((q1_min, q1_max))
    ax.set_ylim((en_min, en_max))

    plot_rez_ellipses(ax)
    ax.legend()
    fig.colorbar(img, ax=ax)
    ax.set_title(
        f"total intensity = {total_intent:.3f}"
        + f"\nCconvolution for {n_measrement} points completed in {t1-t0:.3f} s with {num_worker:1d} cores"
    )

    plt.tight_layout()
    plt.show()
