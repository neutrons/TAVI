from time import time

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse

# from numba import njit, prange


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

    for i in (2, 1, 0):
        if not i == idx:
            mat = quadric_proj(mat, i)

    return 1 / np.sqrt(np.abs(mat[0, 0]))


def model_disp(vq1, vq2):
    """return energy for given Q points
    3d FM J=-1 meV S=1, en=6*S*J*(1-cos(Q))
    """

    sj = 5
    gamma_q = np.cos(2 * np.pi * vq1)
    # gamma_q = (np.cos(2 * np.pi * vq1) + np.cos(2 * np.pi * vq2) + np.cos(2 * np.pi * vq3)) / 3

    disp = 2 * sj * (1 - gamma_q)
    # disp = np.array((disp - 2, disp + 2))

    # reshape if only one band
    num_disp = len(disp.shape)
    if num_disp == 1:
        disp = np.reshape(disp, (1, np.size(disp)))
    return disp


def model_inten(vq1, vq2):
    """return intensity for given Q points
    3d FM J=-1 meV S=1, inten = S/2 for all Qs
    """
    # inten = np.ones_like(vq1, dtype=float) / 2
    # inten = np.array((inten, inten))
    sigma = 0.05
    inten = np.exp(-((vq1 - 0.25) ** 2 + vq2**2) / (2 * sigma**2)) / (2 * np.pi * sigma**2)

    # reshape if only one band
    num_inten = len(inten.shape)
    if num_inten == 1:
        inten = np.reshape(inten, (1, np.size(inten)))

    return inten


def rotation_matrix_3d(theta_deg):
    theta = np.radians(theta_deg)
    c = np.cos(theta)
    s = np.sin(theta)
    mat = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    mat = np.array([[c, 0, -s], [0, 1, 0], [s, 0, c]]) @ mat
    return mat


def resolution_matrix(qx0, qy0, en0):
    """Fake resoltuion matrix mat and prefactor r0
    r0 is a constant, rez_mat is a symmatric positive 4 by 4 matrix
    """

    # sz = np.shape(qx0)

    sigma1, sigma2 = 0.3, 0.02
    sigma3 = 1
    angle = -80
    mat = np.zeros((3, 3))
    mat[0, 0] = 1 / sigma1**2
    mat[1, 1] = 1 / sigma2**2
    mat[2, 2] = 1 / sigma3**2

    rot = rotation_matrix_3d(angle)
    rez_mat = rot.T @ mat @ rot
    r0 = 1
    return r0, rez_mat


def plot_rez_ellipses(ax, xy=(0.25, 0), sigmas=(0.3, 0.02), angle=80, c="k"):
    sigma1, sigma2 = sigmas
    for i in range(3):
        ax.add_artist(
            Ellipse(
                xy=xy,
                width=sigma1 * 2 * (i + 1),
                height=sigma2 * 2 * (i + 1),
                angle=angle,
                edgecolor=c,
                facecolor="none",
                label=f"{i + 1}-sigma",
            )
        )


def generate_pts(num_of_sigmas=3, pts_q=5):
    x = np.linspace(-num_of_sigmas, num_of_sigmas, pts_q + 1)
    v1, v2 = np.meshgrid(x, x, indexing="ij")
    # -------- keep all points ---------
    # vs = np.asarray([np.ravel(v) for v in (v1, v2, v3)])
    # pts_norm = np.asarray(vs)
    # -------- cut the corners based on Gaussian -------
    # g = np.exp(-(v1**2 + v2**2 + v3**2) / 2)
    # idx = g > 1e-4
    # -------- cut the corners based on distance --------
    r_sq = v1**2 + v2**2
    idx = r_sq < num_of_sigmas**2
    # ----------------------------------------
    pts = np.asarray((v1[idx], v2[idx]))
    return pts


def convolution(pts_norm):
    vqh, vqk = trans_mat @ pts_norm
    # ----------------------------------------------------
    # calculate dispersion nergy
    # ----------------------------------------------------
    disp = model_disp(vqh + qh, vqk + qk)
    num_bands, num_pts = disp.shape
    # ----------------------------------------------------
    # Calculate weight from resolution function
    # ----------------------------------------------------
    vq = np.asarray((vqh, vqk))  # shape: (3, num_pts)
    vqe = np.empty((3, num_bands, num_pts))
    vqe[0:2] = vq[:, None, :]
    vqe[2] = disp - en
    prod = np.einsum("ijk,il,ljk->jk", vqe, mat, vqe)
    weights = np.exp(-prod / 2)
    # ----------------------------------------------------
    # trim the corners
    # ----------------------------------------------------
    cut_off = 1e-6
    idx_all = weights > cut_off
    idx = np.any(idx_all, axis=0)
    vq_filtered = vq[:, idx]
    # ----------------------------------------------------
    # calculate intensity
    # ----------------------------------------------------
    elem_vols = elem_vols_init * num_pts_q / np.shape(pts_norm)[1]
    # normalization
    inten = model_inten(vq_filtered[0] + qh, vq_filtered[1] + qk)
    inten_sum = np.sum(inten * weights[:, idx]) * elem_vols
    inten = r0 * inten_sum * np.sqrt(det) / (2 * np.pi) ** (3 / 2)

    return ((vqh, vqk), vq_filtered, inten)


def shift_vector(vec=None):
    "generate the shift vector based on the previous one"
    if vec is None:  # first one
        return np.array((1 / 2, 1 / 2))
    mat = np.array(((-1 / 2, -1 / 2), (1 / 2, -1 / 2)))
    return mat @ vec


if __name__ == "__main__":
    qh, qk, en = (0.25,), (0,), (10,)

    t0 = time()
    # ----------------------------------------------------
    # calculate resolution matrix for all points
    # ----------------------------------------------------
    r0, mat = resolution_matrix(qh, qk, en)
    det = np.linalg.det(mat)
    mat_hkl = quadric_proj(mat, 2)
    # ----------------------------------------------------
    # calculate the incoherent sigmas for E directions
    # ----------------------------------------------------
    sigma_qh = incoh_sigma(mat, 0)
    sigma_qk = incoh_sigma(mat, 1)
    sigma_en = incoh_sigma(mat, 2)
    num_of_sigmas = 3
    min_en, max_en = en - num_of_sigmas * sigma_en, en + num_of_sigmas * sigma_en
    min_qh, max_qh = qh - num_of_sigmas * sigma_qh, qh + num_of_sigmas * sigma_qh
    min_qk, max_qk = qk - num_of_sigmas * sigma_qk, qk + num_of_sigmas * sigma_qk
    # ----------------------------------------------------
    # similarity transformation
    # ----------------------------------------------------
    eigenvalues, eigenvectors = np.linalg.eig(mat_hkl)
    eval_inv_sqrt = 1 / np.sqrt(eigenvalues)
    angle = np.degrees(np.arctan(eigenvectors[0, 1] / eigenvectors[0, 0]))
    # transpose first eignvector!!
    trans_mat = eigenvectors.T @ np.diag(1 / np.sqrt(eigenvalues)) @ eigenvectors
    # ----------------------------------------------------
    # generating points
    # ----------------------------------------------------
    pts_q = 40
    pts_norm = generate_pts(num_of_sigmas, pts_q)
    step_q = 2 * num_of_sigmas / pts_q
    elem_vols_init = step_q**2 * np.prod(eval_inv_sqrt)
    _, num_pts_q = np.shape(pts_norm)

    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(12, 8))

    #  ----------------- plot 1 -----------------
    ax = axes[0, 0]
    c = plt.Circle((0, 0), 3, edgecolor="k", facecolor="none")
    ax.add_artist(c)
    ax.set_xlim((-3, 3))
    ax.set_ylim((-3, 3))
    ax.set_xlabel("Q1")
    ax.set_ylabel("Q2")
    ax.grid(alpha=0.6)
    ax.set_title("initial sampling")

    #  ----------------- plot 2 -----------------
    ax = axes[0, 1]
    plot_rez_ellipses(ax, xy=(0.25, 0), sigmas=eval_inv_sqrt, angle=angle, c="k")
    ax.set_xlim((min_qh, max_qh))
    ax.set_ylim((min_qk, max_qk))
    ax.set_xlabel("Q1")
    ax.set_ylabel("Q2")
    ax.grid(alpha=0.6)
    # plt.colorbar(im)
    ax.set_aspect("equal")
    ax.set_title("similarity transformed")
    #  ----------------- plot 2 -----------------

    ax = axes[0, 2]
    plot_rez_ellipses(ax, xy=(0.25, 0), sigmas=eval_inv_sqrt, angle=angle, c="k")
    ax.set_xlim((min_qh, max_qh))
    ax.set_ylim((min_qk, max_qk))
    ax.set_xlabel("Q1")
    ax.set_ylabel("Q2")
    ax.grid(alpha=0.6)
    ax.set_title("filtered on resolution weight")

    # ----------------- Q1 vs E-----------------
    ax = axes[1, 0]
    ax.plot(x := np.linspace(min_qh, max_qh, 21), model_disp(x, np.ones_like(x) * qk), "-k")
    mat_he = quadric_proj(mat, 1)
    eigenvalues, eigenvectors = np.linalg.eig(mat_he)
    eval_inv_sqrt = 1 / np.sqrt(eigenvalues)
    angle = np.degrees(np.arctan(eigenvectors[0, 1] / eigenvectors[0, 0]))
    plot_rez_ellipses(ax, xy=(qh, model_disp(qh[0], qk[0])), sigmas=eval_inv_sqrt, angle=angle, c="k")
    ax.set_xlim((min_qh, max_qh))
    # sax.set_ylim((min_en, max_en))
    ax.set_ylim((0, 20))
    ax.set_xlabel("Q1")
    ax.set_ylabel("En")
    ax.grid(alpha=0.6)
    ax.set_title("Q1 vs En")

    # ----------------- Q2 vs E-----------------
    ax = axes[1, 1]
    # plot projected resolutio ellipses
    ax.plot(x := np.linspace(min_qk, max_qk, 21), model_disp(np.ones_like(x) * qh, x), "-k")
    mat_ke = quadric_proj(mat, 0)
    eigenvalues, eigenvectors = np.linalg.eig(mat_ke)
    eval_inv_sqrt = 1 / np.sqrt(eigenvalues)
    angle = np.degrees(np.arctan(eigenvectors[0, 1] / eigenvectors[0, 0]))
    plot_rez_ellipses(ax, xy=(qk, model_disp(qh[0], qk[0])), sigmas=eval_inv_sqrt, angle=angle, c="k")
    ax.set_xlim((min_qk, max_qk))
    # ax.set_ylim((min_en, max_en))
    ax.set_ylim((0, 20))
    ax.set_xlabel("Q2")
    ax.set_ylabel("En")
    ax.grid(alpha=0.6)
    ax.set_title("Q2 vs En")
    # ----------------- evolution of intensity-----------------
    ax = axes[1, 2]
    # ax.set_xlim((0, 5))
    # ax.set_ylim((0.032, 0.04))
    ax.set_xlabel("rounds")
    ax.set_ylabel("intensiry")
    ax.grid(alpha=0.6)
    ax.set_title("evolution of intensity")

    # --------------------------------------------------------
    # perform calculation below
    # ---------------------------------------------------------
    (vqh, vqk), vq_filtered, inten = convolution(pts_norm)

    axes[0, 0].plot(pts_norm[0, :], pts_norm[1, :], ".")
    axes[0, 1].plot(vqh + qh, vqk + qk, ".")
    axes[0, 2].plot(vq_filtered[0] + qh, vq_filtered[1] + qk, ".")
    axes[1, 0].plot(vq_filtered[0] + qh, model_disp(vq_filtered[0] + qh, vq_filtered[1] + qk)[0], ".")
    axes[1, 1].plot(vq_filtered[1] + qk, model_disp(vq_filtered[0] + qh, vq_filtered[1] + qk)[0], ".")
    idx_list = [1]
    inten_list = [inten]

    vec = None
    n_round = 1
    sampled_enough = False
    while not sampled_enough:
        vec = shift_vector(vec)
        pts_norm_shifted = pts_norm + (vec * step_q)[:, None]
        (vqh, vqk), vq_filtered, inten_new = convolution(pts_norm_shifted)

        axes[0, 0].plot(pts_norm_shifted[0, :], pts_norm_shifted[1, :], ".")
        axes[0, 1].plot(vqh + qh, vqk + qk, ".")
        axes[0, 2].plot(vq_filtered[0] + qh, vq_filtered[1] + qk, ".")
        axes[1, 0].plot(vq_filtered[0] + qh, model_disp(vq_filtered[0] + qh, vq_filtered[1] + qk)[0], ".")
        axes[1, 1].plot(vq_filtered[1] + qk, model_disp(vq_filtered[0] + qh, vq_filtered[1] + qk)[0], ".")
        idx_list.append(n_round + 1)
        if np.abs(inten_new - inten) / inten < 1e-2:
            sampled_enough = True
        else:
            pts_norm = np.concatenate((pts_norm, pts_norm_shifted), axis=1)
            n_round += 1
        inten = (inten + inten_new) / 2
        inten_list.append(inten)

    # pts_norm_2 = pts_norm + np.asarray((step_q / 2, step_q / 2))[:, None]
    # pts_norm_3 = np.concatenate((pts_norm, pts_norm_2), axis=1) + np.asarray((-step_q / 2, 0))[:, None]
    # pts_norm_4 = (
    #     np.concatenate((pts_norm, pts_norm_2, pts_norm_3), axis=1) + np.asarray((-step_q / 4, -step_q / 4))[:, None]
    # )

    axes[1, 2].plot(idx_list, inten_list, "-o")

    plt.suptitle("Adaptive Sampling Demo (3D)")
    plt.tight_layout()
    plt.show()
