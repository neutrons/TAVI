import concurrent.futures
from time import time

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse
from numba import njit, prange


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
            mat = quadric_proj(mat, i)

    return 1 / np.sqrt(np.abs(mat[0, 0]))


# @njit(parallel=True)
def model_disp(vq1, vq2, vq3):
    """return energy for given Q points
    3d FM J=-1 meV S=1, en=6*S*J*(1-cos(Q))
    """

    sj = 5
    gamma_q = np.cos(2 * np.pi * vq1)
    # gamma_q = (np.cos(2 * np.pi * vq1) + np.cos(2 * np.pi * vq2) + np.cos(2 * np.pi * vq3)) / 3

    disp = 2 * sj * (1 - gamma_q)
    disp = np.array((disp - 2, disp + 2))

    # reshape if only one band
    num_disp = len(disp.shape)
    if num_disp == 1:
        disp = np.reshape(disp, (1, np.size(disp)))
    return disp


# @njit(parallel=True)
def model_inten(vq1, vq2, vq3):
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


def resolution_matrix(qx0, qy0, qz0, en0):
    """Fake resoltuion matrix mat and prefactor r0
    r0 is a constant, rez_mat is a symmatric positive 4 by 4 matrix
    """
    # sz = np.shape(qx0)

    sigma1, sigma2 = 0.3, 0.02
    sigma3 = sigma4 = 1
    angle = -80
    mat = np.zeros((4, 4))
    mat[0, 0] = 1 / sigma1**2
    mat[1, 1] = 1 / sigma3**2
    mat[2, 2] = 1 / sigma4**2
    mat[3, 3] = 1 / sigma2**2

    rot = rotation_matrix_4d(angle)
    rez_mat = rot.T @ mat @ rot
    r0 = 1
    return r0, rez_mat


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


@njit(parallel=True, nogil=True)
def compute_weights(vqe, mat):
    _, num_bands, num_pts = vqe.shape
    weights = np.empty((num_bands, num_pts))

    for i in prange(num_bands):
        for j in range(num_pts):
            v = vqe[:, i, j]  # shape: (4,)
            tmp = 0.0
            for k in range(4):
                for l in range(4):
                    tmp += v[k] * mat[k, l] * v[l]
            weights[i, j] = np.exp(-0.5 * tmp)

    return weights


def convolution(qh, qk, ql, en):
    # ----------------------------------------------------
    # calculate resolution matrix for all points
    # ----------------------------------------------------
    r0, mat = resolution_matrix(qh, qk, ql, en)
    mat_hkl = quadric_proj(mat, 3)
    # ----------------------------------------------------
    # calculate the incoherent sigmas for E directions
    # ----------------------------------------------------
    sigma_en = incoh_sigma(mat, 3)
    num_of_sigmas = 3
    min_en, max_en = en - num_of_sigmas * sigma_en, en + num_of_sigmas * sigma_en
    # ----------------------------------------------------
    # similarity transformation
    # ----------------------------------------------------
    eigenvalues, eigenvectors = np.linalg.eig(mat_hkl)
    eval_inv_sqrt = 1 / np.sqrt(eigenvalues)
    trans_mat = eigenvectors * eval_inv_sqrt * eigenvectors
    # ----------------------------------------------------
    # start sampling
    # ----------------------------------------------------
    pts_q = 5
    sampled_enough = False
    while not sampled_enough:
        x = np.linspace(-num_of_sigmas, num_of_sigmas, pts_q + 1)
        elem_vols = (2 * num_of_sigmas / pts_q) ** 3
        v1, v2, v3 = np.meshgrid(x, x, x, indexing="ij")
        vs = np.asarray([np.ravel(v) for v in (v1, v2, v3)])
        pts_norm = np.asarray(vs)

        # cut the corners beyond 3*sigma when dimemsion is higher than 1
        # g = np.exp(-(v1**2 + v2**2 + v3**2) / 2)
        # idx = g > 1e-4
        # pts_norm = np.asarray((v1[idx], v2[idx], v3[idx]))

        vqh, vqk, vql = trans_mat @ pts_norm
        vqh += qh
        vqk += qk
        vql += ql
        elem_vols *= np.prod(eval_inv_sqrt)
        # ----------------------------------------------------
        # determine if sampled enough based on steps along energy
        # ----------------------------------------------------
        disp = model_disp(vqh, vqk, vql)
        num_bands, num_pts = disp.shape
        # ----------------------------------------------------
        # get rid of the ones that are not in the ellipsoid
        # ----------------------------------------------------
        max_disp, min_disp = np.max(disp), np.min(disp)
        if max_disp < min_en or min_disp > max_en:
            return 0.0  # zero intensity
        # ----------------------------------------------------
        # determine if sampled enough based on steps along energy
        # ----------------------------------------------------
        max_en_step = (max_disp - min_disp) / pts_q * 2
        ratio = max_en_step / sigma_en / 2
        if ratio > 1:
            pts_q = int(pts_q * ratio) + 1
            print(f"ratio={ratio:.2f} for (Q,E)=({qh:.2f}, {en:.2f})")
            continue
        break

    # ----------------------------------------------------
    # Enough sampled. Calculate weight from resolution function
    # ----------------------------------------------------
    vq = np.asarray((vqh - qh, vqk - qk, vql - ql))  # shape: (3, num_pts)
    vqe = np.empty((4, num_bands, num_pts))
    vqe[0:3] = vq[:, None, :]
    vqe[3] = disp - en
    # prod = np.einsum("ijk,il,ljk->jk", vqe, mat, vqe)
    # weights = np.exp(-prod / 2)  # shape: (num_bands, num_pts)
    # ------------------------------------
    weights = compute_weights(vqe, mat)
    # ------------------------------------

    # don't bother if the weight is already too small
    if np.max(weights) < 1e-6:
        return 0.0  # zero intensity

    # ----------------------------------------------------
    # trim the corners
    # ----------------------------------------------------
    cut_off = 1e-6
    idx_all = weights > cut_off
    idx = np.any(idx_all, axis=0)
    # percent = (np.size(idx) - np.count_nonzero(idx)) / np.size(idx) * 100
    # print(f"{percent:.2f}% of points discarded.")

    # all small weights because dispersion parallel to ellipsoid
    if not np.any(idx_all):
        return 0.0  # zero intensity

    vq_filtered = vq[:, idx]
    weights_filtered = weights[:, idx]
    inten = model_inten(*vq_filtered)

    # normalization
    det = np.linalg.det(mat)
    inten_sum = np.sum(inten * weights_filtered) * elem_vols
    return r0 * inten_sum * np.sqrt(det) / (2 * np.pi) ** 2


# @njit(parallel=True)
# def parallel_convolution(qe_mesh):
#     n_pts = qe_mesh.shape[1]
#     results = np.empty(n_pts)
#     for i in prange(n_pts):
#         qh, qk, ql, en = qe_mesh[0, i], qe_mesh[1, i], qe_mesh[2, i], qe_mesh[3, i]
#         results[i] = convolution(qh, qk, ql, en)
#     return results


if __name__ == "__main__":
    # ----------------------------------------------------
    # points being measured
    # qe_mesh has the dimension (4, n_pts_of_measurement)
    # flatten for meshed measurement
    # ----------------------------------------------------
    q1_min, q1_max, q1_step = -1, 1, 0.02
    en_min, en_max, en_step = -3, 25, 0.2
    q2 = 0
    q3 = 0

    q1 = np.linspace(q1_min, q1_max, int((q1_max - q1_min) / q1_step) + 1)
    en = np.linspace(en_min, en_max, int((en_max - en_min) / en_step) + 1)
    vq1, vq2, vq3, ven = np.meshgrid(q1, q2, q3, en, indexing="ij")

    sz = np.shape(vq1)  # sz = (n_q1, n_q2, n_q3 , n_en)
    qe_mesh = np.asarray([np.ravel(v) for v in (vq1, vq2, vq3, ven)])

    t0 = time()
    num_worker = 4
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_worker) as executor:
        results = executor.map(convolution, *qe_mesh)
    measurement_inten = np.asarray(list(results)).reshape(sz)

    # qe_mesh = np.ascontiguousarray(qe_mesh)
    # measurement_inten = parallel_convolution(qe_mesh).reshape(sz)
    print(f"Convolution completed in {(t1:=time())-t0:.4f} s")
    # total intensity should be close to S/2 *(q1_max - q1_min) * 2p*i
    total_intent = np.sum(measurement_inten) * q1_step * en_step / (q1_max - q1_min)

    # ----------------------------------------------------
    # plot 2D contour
    # ----------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 4))
    idx = np.s_[:, 0, 0, :]
    img = ax.pcolormesh(vq1[idx], ven[idx], measurement_inten[idx], cmap="turbo", vmin=0, vmax=0.5)

    ax.grid(alpha=0.6)
    ax.set_xlabel("Q1")
    ax.set_ylabel("En")
    ax.set_xlim((q1_min, q1_max))
    ax.set_ylim((en_min, en_max))

    plot_rez_ellipses(ax)
    disp = model_disp(q1, np.zeros_like(q1), np.zeros_like(q1))
    for i in range(np.shape(disp)[0]):
        ax.plot(q1, disp[i], "-w")

    ax.legend()
    fig.colorbar(img, ax=ax)
    ax.set_title(
        f"1D FM chain S=1 J=-5, total intensity = {total_intent:.3f}"
        + f"\n3D Convolution for {np.shape(qe_mesh)[1]} points completed in {t1-t0:.3f} s with {num_worker:1d} cores"
    )

    plt.tight_layout()
    plt.show()
