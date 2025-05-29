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


def coh_sigma(mat, axis):
    """Coherent sigma"""
    idx = int(axis)

    return 1 / np.sqrt(np.abs(mat[idx, idx]))


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
    sigma3 = sigma4 = 0.2
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
                label=f"{i + 1}-sigma",
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
            weights[i, j] = tmp
            # weights[i, j] = np.exp(-0.5 * tmp)

    return weights


def generate_pts(sigma_qs, num_of_sigmas, num_pts, mat_hkl):
    (sigma_qh_incoh, sigma_qk_incoh, sigma_ql_incoh) = sigma_qs
    sigma_qh_range, sigma_qk_range, sigma_ql_range = (
        num_of_sigmas * sigma_qh_incoh,
        num_of_sigmas * sigma_qk_incoh,
        num_of_sigmas * sigma_ql_incoh,
    )
    pts_qh, pts_qk, pts_ql = num_pts
    qh_list = np.linspace(-sigma_qh_range, sigma_qh_range, pts_qh + 1)
    qk_list = np.linspace(-sigma_qk_range, sigma_qk_range, pts_qk + 1)
    ql_list = np.linspace(-sigma_ql_range, sigma_ql_range, pts_ql + 1)
    vq = np.meshgrid(qh_list, qk_list, ql_list, indexing="ij")
    # -------- cut the corners based on distance --------
    r_sq = np.einsum("ijkl,im,mjkl->jkl", vq, mat_hkl, vq)
    idx = r_sq < num_of_sigmas**2

    return vq, idx


def get_max_step(arr, axis: int):
    diff_arr = np.abs(np.diff(arr, axis=axis))
    # diff_idx = np.all(np.isnan(diff_arr), axis=axis)
    step = np.nanmean(diff_arr, axis=axis)
    # print(f"max_step={step:.3f}")
    return np.nanmax(step)


def convolution(qh, qk, ql, en):
    # ----------------------------------------------------
    # calculate resolution matrix for all points
    # ----------------------------------------------------
    r0, mat = resolution_matrix(qh, qk, ql, en)
    mat_hkl = quadric_proj(mat, 3)
    # ----------------------------------------------------
    # calculate the incoherent sigmas for all Q and E directions
    # ----------------------------------------------------
    sigma_qs = tuple(incoh_sigma(mat, i) for i in range(3))
    sigma_en_incoh = incoh_sigma(mat, 3)
    num_of_sigmas = 3
    min_en, max_en = en - num_of_sigmas * sigma_en_incoh, en + num_of_sigmas * sigma_en_incoh
    sigma_en_coh = coh_sigma(mat, 3)
    en_rez = sigma_en_coh / 5

    # ----------------------------------------------------
    # Calculate elemental volume
    # ----------------------------------------------------
    eigenvalues = np.linalg.eigvalsh(mat_hkl)
    eval_inv_sqrt = 1 / np.sqrt(eigenvalues)
    elem_vols = np.prod(eval_inv_sqrt) * (2 * num_of_sigmas) ** 3
    # ----------------------------------------------------
    # Adaptive sampling
    # ----------------------------------------------------
    pts = [10, 10, 10]
    sampled_enough = False

    while not sampled_enough:
        (vqh, vqk, vql), idx = generate_pts(sigma_qs, num_of_sigmas, pts, mat_hkl)
        # ----------------------------------------------------
        # determine if sampled enough based on steps along energy
        # ----------------------------------------------------
        disp = model_disp(vqh[idx] + qh, vqk[idx] + qk, vql[idx] + ql)
        num_bands, num_pts = disp.shape

        # Skip if all dispersion is outside the relevant energy window
        if np.max(disp) < min_en or np.min(disp) > max_en:
            return 0.0
        # ----------------------------------------------------
        # determine if sampled enough based on steps along energy
        # ----------------------------------------------------
        disp_arr = np.full(shape=(num_bands,) + idx.shape, fill_value=np.nan)
        disp_arr[(slice(None),) + np.nonzero(idx)] = disp

        # Compute max energy steps
        steps = [get_max_step(disp_arr, axis=i) for i in (1, 2, 3)]

        # Adaptive sampling
        sampled_enough_flags = []
        for i, (step, pt) in enumerate(zip(steps, pts)):
            if step > en_rez:
                pts[i] = pt * (int(step / en_rez) + 1)
                sampled_enough_flags.append(False)
            else:
                sampled_enough_flags.append(True)

        sampled_enough = all(sampled_enough_flags)

        if sampled_enough:
            print(f"Enough sample for (Q1, E) = ({qh:.2f}, {en:.2f}), num of pts = {pts}")
            print(f"steps in energy = ({steps[0]:.2f}, {steps[1]:.2f}, {steps[2]:.2f})")
            vqh, vqk, vql = vqh[idx], vqk[idx], vql[idx]
            break

    # ----------------------------------------------------
    # Enough sampled. Calculate weight from resolution function
    # ----------------------------------------------------
    elem_vols /= np.prod(pts)

    vq = np.array((vqh, vqk, vql))  # shape: (3, num_pts)
    vqe = np.empty((4, num_bands, num_pts))
    vqe[0:3] = vq[:, None, :]
    vqe[3] = disp - en

    weights = compute_weights(vqe, mat)  # shape: (num_bands, num_pts)
    # ------------------------------------

    # # don't bother if the weight is already too small
    # # if np.max(weights) < 1e-6:
    # #     return 0.0  # zero intensity

    # if np.min(weights) > 9:  # out of 3 sigma squared
    #     return 0.0  # zero intensity

    # # ----------------------------------------------------
    # # trim the corners
    # # ----------------------------------------------------
    # # TODO
    # cut_off = 9  # 3 sigma squared
    # idx_all = weights < cut_off
    # idx = np.any(idx_all, axis=0)
    # # percent = (np.size(idx) - np.count_nonzero(idx)) / np.size(idx) * 100
    # # print(f"{percent:.2f}% of points discarded.")

    # # all small weights because dispersion parallel to ellipsoid
    # if not np.any(idx_all):
    #     return 0.0  # zero intensity

    # Fast exit if everything is outside the 3-sigma volume
    if np.min(weights) > 9 or not np.any(weights < 9):
        return 0.0

    idx_keep = np.any(weights < 9, axis=0)
    vq_filtered = vq[:, idx_keep]
    weights_filtered = np.exp(-weights[:, idx_keep] / 2)
    inten = model_inten(*vq_filtered)

    # normalization
    det = np.linalg.det(mat)
    inten_sum = np.sum(inten * weights_filtered) * elem_vols
    return r0 * inten_sum * np.sqrt(det) / (2 * np.pi) ** 2


if __name__ == "__main__":
    # ----------------------------------------------------
    # points being measured
    # qe_mesh has the dimension (4, n_pts_of_measurement)
    # flatten for meshed measurement
    # ----------------------------------------------------
    q1_min, q1_max, q1_step = 0, 1, 0.02
    en_min, en_max, en_step = -3, 25, 0.2
    q2 = 0
    q3 = 0

    q1 = np.linspace(q1_min, q1_max, int((q1_max - q1_min) / q1_step) + 1)
    en = np.linspace(en_min, en_max, int((en_max - en_min) / en_step) + 1)
    vq1, vq2, vq3, ven = np.meshgrid(q1, q2, q3, en, indexing="ij")

    sz = np.shape(vq1)  # sz = (n_q1, n_q2, n_q3 , n_en)
    qe_mesh = np.asarray([np.ravel(v) for v in (vq1, vq2, vq3, ven)])

    t0 = time()
    num_worker = 8
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_worker) as executor:
        results = executor.map(convolution, *qe_mesh)
    measurement_inten = np.asarray(list(results)).reshape(sz)

    # qe_mesh = np.ascontiguousarray(qe_mesh)
    # measurement_inten = parallel_convolution(qe_mesh).reshape(sz)
    print(f"Convolution completed in {(t1 := time()) - t0:.4f} s")
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
        + f"\n3D Convolution for {np.shape(qe_mesh)[1]} points completed in {t1 - t0:.3f} s with {num_worker:1d} cores"
    )

    plt.tight_layout()
    plt.show()
