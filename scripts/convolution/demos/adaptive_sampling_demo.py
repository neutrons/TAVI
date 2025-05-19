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
    inten = np.ones_like(vq1, dtype=float) / 2
    inten = np.array((inten, inten))

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
                label=f"{i+1}-sigma",
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
    pts_norm = np.asarray((v1[idx], v2[idx]))
    return pts_norm


# @njit(parallel=True, nogil=True)
# def compute_weights(vqe, mat):
#     _, num_bands, num_pts = vqe.shape
#     weights = np.empty((num_bands, num_pts))

#     for i in prange(num_bands):
#         for j in range(num_pts):
#             v = vqe[:, i, j]  # shape: (4,)
#             tmp = 0.0
#             for k in range(4):
#                 for l in range(4):
#                     tmp += v[k] * mat[k, l] * v[l]
#             weights[i, j] = np.exp(-0.5 * tmp)

#     return weights


if __name__ == "__main__":
    qh, qk, en = (0.25,), (0,), (10,)

    t0 = time()
    # ----------------------------------------------------
    # calculate resolution matrix for all points
    # ----------------------------------------------------
    r0, mat = resolution_matrix(qh, qk, en)
    det = np.linalg.det(mat)
    mat_hkl = quadric_proj(mat, -1)
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
    # generating points
    # ----------------------------------------------------
    pts_q = 10
    pts_norm = generate_pts(num_of_sigmas, pts_q)
    step_q = 2 * num_of_sigmas / pts_q
    elem_vols = step_q**2
    pts_norm_2 = pts_norm + np.asarray((step_q / 2, step_q / 2))[:, None]
    pts_norm_3 = np.concatenate((pts_norm, pts_norm_2), axis=1) + np.asarray((-step_q / 2, 0))[:, None]
    pts_norm_4 = (
        np.concatenate((pts_norm, pts_norm_2, pts_norm_3), axis=1) + np.asarray((-step_q / 4, -step_q / 4))[:, None]
    )
    # ----------------------------------------------------
    # similarity transformation
    # ----------------------------------------------------
    eigenvalues, eigenvectors = np.linalg.eig(mat_hkl)
    eval_inv_sqrt = 1 / np.sqrt(eigenvalues)
    elem_vols *= np.prod(eval_inv_sqrt)
    angle = np.degrees(np.arctan(eigenvectors[0, 1] / eigenvectors[0, 0]))
    # transpose first eignvector!!
    trans_mat = eigenvectors.T @ np.diag(1 / np.sqrt(eigenvalues)) @ eigenvectors

    vqh, vqk = trans_mat @ pts_norm
    vqh_2, vqk_2 = trans_mat @ pts_norm_2
    vqh_3, vqk_3 = trans_mat @ pts_norm_3
    vqh_4, vqk_4 = trans_mat @ pts_norm_4
    # ----------------------------------------------------
    # determine if sampled enough based on steps along energy
    # ----------------------------------------------------
    disp = model_disp(vqh + qh, vqk + qk)
    num_bands, num_pts = disp.shape
    disp_2 = model_disp(vqh_2 + qh, vqk_2 + qk)
    disp_3 = model_disp(vqh_3 + qh, vqk_3 + qk)
    disp_4 = model_disp(vqh_4 + qh, vqk_4 + qk)
    # ----------------------------------------------------
    # Calculate weight from resolution function
    # ----------------------------------------------------
    vq = np.asarray((vqh, vqk))  # shape: (3, num_pts)
    vqe = np.empty((3, num_bands, num_pts))
    vqe[0:2] = vq[:, None, :]
    vqe[2] = disp - en
    prod = np.einsum("ijk,il,ljk->jk", vqe, mat, vqe)
    weights = np.exp(-prod / 2)
    # ------------------------------------
    # weights = compute_weights(vqe, mat)  # shape: (num_bands, num_pts)

    # --------------
    vq_2 = np.asarray((vqh_2, vqk_2))  # shape: (3, num_pts)
    vqe_2 = np.empty((3, num_bands, num_pts))
    vqe_2[0:2] = vq_2[:, None, :]
    vqe_2[2] = disp_2 - en
    # weights_2 = compute_weights(vqe_2, mat)  # shape: (num_bands, num_pts)
    prod = np.einsum("ijk,il,ljk->jk", vqe_2, mat, vqe_2)
    weights_2 = np.exp(-prod / 2)
    # --------------
    vq_3 = np.asarray((vqh_3, vqk_3))  # shape: (3, num_pts)
    vqe_3 = np.empty((3, num_bands, num_pts * 2))
    vqe_3[0:2] = vq_3[:, None, :]
    vqe_3[2] = disp_3 - en
    # weights_3 = compute_weights(vqe_3, mat)  # shape: (num_bands, num_pts)
    prod = np.einsum("ijk,il,ljk->jk", vqe_3, mat, vqe_3)
    weights_3 = np.exp(-prod / 2)
    # --------------
    vq_4 = np.asarray((vqh_4, vqk_4))  # shape: (3, num_pts)
    vqe_4 = np.empty((3, num_bands, num_pts * 4))
    vqe_4[0:2] = vq_4[:, None, :]
    vqe_4[2] = disp_4 - en
    # weights_4 = compute_weights(vqe_4, mat)  # shape: (num_bands, num_pts)
    prod = np.einsum("ijk,il,ljk->jk", vqe_4, mat, vqe_4)
    weights_4 = np.exp(-prod / 2)

    # ----------------------------------------------------
    # trim the corners
    # ----------------------------------------------------
    cut_off = 1e-9
    idx_all = weights > cut_off
    idx = np.any(idx_all, axis=0)
    vq_filtered = vq[:, idx]
    # weights_filtered = weights[:, idx]
    idx_all_2 = weights_2 > cut_off
    idx_2 = np.any(idx_all_2, axis=0)
    vq_filtered_2 = vq_2[:, idx_2]
    # -----------
    idx_all_3 = weights_3 > cut_off
    idx_3 = np.any(idx_all_3, axis=0)
    vq_filtered_3 = vq_3[:, idx_3]
    # -----------
    idx_all_4 = weights_4 > cut_off
    idx_4 = np.any(idx_all_4, axis=0)
    vq_filtered_4 = vq_4[:, idx_4]

    # ----------------------------------------------------
    # calculate intensity
    # ----------------------------------------------------
    elem_vols_2 = elem_vols
    elem_vols_3 = elem_vols / 2
    elem_vols_4 = elem_vols / 4
    # normalization
    inten = model_inten(*vq_filtered)
    inten_sum = np.sum(inten * weights[:, idx]) * elem_vols
    inten = r0 * inten_sum * np.sqrt(det) / (2 * np.pi) ** (3 / 2)

    inten_2 = model_inten(*vq_filtered_2)
    inten_sum_2 = np.sum(inten_2 * weights_2[:, idx_2]) * elem_vols_2
    inten_2 = r0 * inten_sum_2 * np.sqrt(det) / (2 * np.pi) ** (3 / 2)
    inten_2 = (inten + inten_2) / 2

    inten_3 = model_inten(*vq_filtered_3)
    inten_sum_3 = np.sum(inten_3 * weights_3[:, idx_3]) * elem_vols_3
    inten_3 = r0 * inten_sum_3 * np.sqrt(det) / (2 * np.pi) ** (3 / 2)
    inten_3 = (inten_2 + inten_3) / 2

    inten_4 = model_inten(*vq_filtered_4)
    inten_sum_4 = np.sum(inten_4 * weights_4[:, idx_4]) * elem_vols_4
    inten_4 = r0 * inten_sum_4 * np.sqrt(det) / (2 * np.pi) ** (3 / 2)
    inten_4 = (inten_3 + inten_4) / 2

    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(12, 8))

    #  ----------------- plot 1 -----------------
    ax = axes[0, 0]
    ax.plot(pts_norm[0, :], pts_norm[1, :], ".")
    ax.plot(pts_norm_2[0, :], pts_norm_2[1, :], ".")
    ax.plot(pts_norm_3[0, :], pts_norm_3[1, :], ".")
    ax.plot(pts_norm_4[0, :], pts_norm_4[1, :], ".")
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
    ax.plot(vqh + qh, vqk + qk, ".")
    ax.plot(vqh_2 + qh, vqk_2 + qk, ".")
    ax.plot(vqh_3 + qh, vqk_3 + qk, ".")
    ax.plot(vqh_4 + qh, vqk_4 + qk, ".")
    # make contour of dispersion
    # qh_list = np.linspace(min_qh, max_qh, 201)
    # qk_list = np.linspace(min_qk, max_qk, 201)
    # x, y = np.meshgrid(qh_list, qk_list)
    # im = ax.pcolormesh(x, y, model_disp(x, y), cmap="turbo", vmin=0, vmax=20)
    plot_rez_ellipses(ax, xy=(0.25, 0), sigmas=eval_inv_sqrt, angle=angle, c="k")
    ax.set_xlim((min_qh, max_qh))
    ax.set_ylim((min_qk, max_qk))
    ax.set_xlabel("Q1")
    ax.set_ylabel("Q2")
    ax.grid(alpha=0.6)
    # plt.colorbar(im)
    # ax.set_aspect("equal")
    ax.set_title("similarity transformed")
    #  ----------------- plot 2 -----------------
    ax = axes[0, 2]
    ax.plot(vq_filtered[0] + qh, vq_filtered[1] + qk, ".")
    ax.plot(vq_filtered_2[0] + qh, vq_filtered_2[1] + qk, ".")
    ax.plot(vq_filtered_3[0] + qh, vq_filtered_3[1] + qk, ".")
    ax.plot(vq_filtered_4[0] + qh, vq_filtered_4[1] + qk, ".")

    plot_rez_ellipses(ax, xy=(0.25, 0), sigmas=eval_inv_sqrt, angle=angle, c="k")
    ax.set_xlim((min_qh, max_qh))
    ax.set_ylim((min_qk, max_qk))
    ax.set_xlabel("Q1")
    ax.set_ylabel("Q2")
    ax.grid(alpha=0.6)
    ax.set_title("filtered on resolution weight")

    # ----------------- Q1 vs E-----------------
    ax = axes[1, 0]
    ax.plot(vq_filtered[0] + qh, model_disp(vq_filtered[0] + qh, vq_filtered[1] + qk)[0], "o")
    ax.plot(vq_filtered_2[0] + qh, model_disp(vq_filtered_2[0] + qh, vq_filtered_2[1] + qk)[0], "o")
    ax.plot(vq_filtered_3[0] + qh, model_disp(vq_filtered_3[0] + qh, vq_filtered_3[1] + qk)[0], "o")
    ax.plot(vq_filtered_4[0] + qh, model_disp(vq_filtered_4[0] + qh, vq_filtered_4[1] + qk)[0], "o")
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
    ax.plot(vq_filtered[1] + qk, model_disp(vq_filtered[0] + qh, vq_filtered[1] + qk)[0], "o")
    ax.plot(vq_filtered_2[1] + qk, model_disp(vq_filtered_2[0] + qh, vq_filtered_2[1] + qk)[0], "o")
    ax.plot(vq_filtered_3[1] + qk, model_disp(vq_filtered_3[0] + qh, vq_filtered_3[1] + qk)[0], "o")
    ax.plot(vq_filtered_4[1] + qk, model_disp(vq_filtered_4[0] + qh, vq_filtered_4[1] + qk)[0], "o")
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
    # ----------------- Q2 vs E-----------------
    ax = axes[1, 2]
    ax.plot([1, 2, 3, 4], [inten, inten_2, inten_3, inten_4], "-o")
    ax.set_xlim((0, 5))
    # ax.set_ylim((0.032, 0.04))
    ax.set_xlabel("rounds")
    ax.set_ylabel("intensiry")
    ax.grid(alpha=0.6)
    ax.set_title("evolution of intensity")

    plt.suptitle("Adaptive Sampling Demo (3D)")
    plt.tight_layout()
    plt.show()
