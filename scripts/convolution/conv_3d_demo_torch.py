import functools
from time import time

import matplotlib.pyplot as plt
import numpy as np
import torch
from matplotlib.patches import Ellipse


# -------------------------------------------------------
# user input model_disp and model_inten
# -------------------------------------------------------
def model_disp(vq1, vq2, vq3):
    """return energy for given Q points
    3d FM J=-1 meV S=1, en=6*S*J*(1-cos(Q))
    """

    sj = 5
    # gamma_q = np.cos(2 * np.pi * vq1)
    gamma_q = (np.cos(2 * np.pi * vq1) + np.cos(2 * np.pi * vq2) + np.cos(2 * np.pi * vq3)) / 3

    disp = 2 * sj * (1 - gamma_q)
    disp = np.array((disp - 2, disp + 2))

    # reshape if only one band
    num_disp = len(disp.shape)
    if num_disp == 1:
        disp = np.reshape(disp, (1, np.size(disp)))
    return disp


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


# -------------------------------------------------------
# fake resolution matrix and resolution ellipses
# -------------------------------------------------------
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


def resolution_matrix(hkl, en):
    """Fake resoltuion matrix mat and prefactor r0
    r0 is a constant, rez_mat is a symmetric positive 4 by 4 matrix
    """

    sigma1, sigma2 = 0.3, 0.02
    sigma3 = sigma4 = 0.2
    angle = -80
    mat = np.diag([1 / sigma1**2, 1 / sigma3**2, 1 / sigma4**2, 1 / sigma2**2])

    rot = rotation_matrix_4d(angle)
    rez_mat = rot.T @ mat @ rot
    r0 = 1

    return tuple((hkl[i], en[j], r0, rez_mat) for i in range(np.shape(hkl)[0]) for j in range(np.size(en)))


def plot_rez_ellipses(ax):
    sigma1, sigma2 = 0.3, 0.02
    angle = 80
    for i in range(3):
        ax.add_artist(
            Ellipse(
                xy=(2, 0),
                width=sigma1 * 2 * (i + 1),
                height=sigma2 * 2 * (i + 1),
                angle=angle,
                edgecolor="w",
                facecolor="none",
                label=f"{i + 1}-sigma",
            )
        )


# -------------------------------------------------------
# functions required for resolution convolution
# -------------------------------------------------------


def quadric_proj(quadric: np.ndarray, idx: int) -> np.ndarray:
    """projects along one axis of the quadric

    dimensino of input arry is n by n
    dimensiont of output array is (n-1) by (n-1)
    """

    # delete if orthogonal
    if np.abs(qii := quadric[idx, idx]) < 1e-8:
        mask = np.arange(quadric.shape[0]) != idx
        return quadric[np.ix_(mask, mask)]

    # row/column along which to perform the orthogonal projection
    # symmetrise if not symmetric, normalise to indexed component
    vec = 0.5 * (quadric[idx, :] + quadric[:, idx]) / np.sqrt(qii)
    ortho_proj = quadric - np.outer(vec, vec)  # projected quadric

    # return np.delete(np.delete(ortho_proj, idx, axis=0), idx, axis=1)
    mask = np.arange(ortho_proj.shape[0]) != idx
    return ortho_proj[np.ix_(mask, mask)]


def incoh_sigma_en(mat: np.ndarray) -> float:
    """Incoherent sigma for energy"""

    elem = quadric_proj(quadric_proj(quadric_proj(mat, 2), 1), 0)[0, 0]

    return 1 / np.sqrt(np.abs(elem))


def incoh_sigma_qs(mat: np.ndarray) -> tuple[float, float, float]:
    """Incoherent sigmas for q1, q2 and q3"""

    mat_2 = quadric_proj(mat, 2)
    elem1 = abs(quadric_proj(mat_2, 1)[0, 0])
    elem2 = abs(quadric_proj(mat_2, 0)[0, 0])
    elem3 = abs(quadric_proj(quadric_proj(mat, 1), 0)[0, 0])

    return (1 / np.sqrt(elem1), 1 / np.sqrt(elem2), 1 / np.sqrt(elem3))


def coh_sigma(mat: np.ndarray, axis: int):
    """Coherent sigma along a given axis"""
    idx = int(axis)

    return 1 / np.sqrt(np.abs(mat[idx, idx]))


# @njit(parallel=True, nogil=True)
# def compute_weights(vqe: np.ndarray, mat: np.ndarray) -> np.ndarray:
#     """calculate weiget
#     vqe has shape (4, num_bands, num_pts)
#     mat has shape (4, 4)
#     weights = np.einsum("ijk,il,ljk->jk", vqe, mat_qe, vqe)
#     """
#     _, num_bands, num_pts = vqe.shape
#     weights = np.empty((num_bands, num_pts))

#     for i in prange(num_bands):
#         for j in range(num_pts):
#             v = vqe[:, i, j]  # shape: (4,)
#             tmp = 0.0
#             for k in range(4):
#                 for l in range(4):
#                     tmp += v[k] * mat[k, l] * v[l]
#             weights[i, j] = tmp

#     return weights


def compute_weights_torch(vqe: torch.Tensor, mat: torch.Tensor) -> torch.Tensor:
    """
    vqe: shape (4, num_bands, num_pts)
    mat: shape (4, 4)
    return: shape (num_bands, num_pts)
    """
    # shape: (4, num_bands, num_pts) → (num_bands, num_pts, 4)
    v = vqe.permute(1, 2, 0)

    mv = torch.matmul(v, mat)  # shape: (num_bands, num_pts, 4)
    weights = torch.sum(v * mv, dim=-1)  # shape: (num_bands, num_pts)

    return weights


@functools.cache
def generate_meshgrid(num_of_sigmas=3, num_pts=(10, 10, 10)):
    pts_qh, pts_qk, pts_ql = num_pts
    qh = np.linspace(-num_of_sigmas, num_of_sigmas, pts_qh + 1)
    qk = np.linspace(-num_of_sigmas, num_of_sigmas, pts_qk + 1)
    ql = np.linspace(-num_of_sigmas, num_of_sigmas, pts_ql + 1)
    return np.meshgrid(qh, qk, ql, indexing="ij")  # shape (3, N1, N2, N3)


def generate_pts(sigma_qs, mat_hkl, num_of_sigmas=3, num_pts=(10, 10, 10)):
    """Generate points in a 3D mesh, cut the points at the corners"""
    (sigma_qh_incoh, sigma_qk_incoh, sigma_ql_incoh) = sigma_qs

    vq_h, vq_k, vq_l = generate_meshgrid(num_of_sigmas, num_pts)
    vq = (vq_h * sigma_qh_incoh, vq_k * sigma_qk_incoh, vq_l * sigma_ql_incoh)

    # -------- cut the corners based on distance --------
    r_sq = np.einsum("i...,ij,j...->...", vq, mat_hkl, vq)
    idx = r_sq < num_of_sigmas**2  # Ellipsoid mask
    return (vq[0][idx], vq[1][idx], vq[2][idx]), idx


def get_max_step(arr, axis: int):
    """Get max step along a given axis. Return zero if all NaN"""
    # shape of arr (num_bands, N1, N2, N3)
    diff_arr = np.abs(np.diff(arr, axis=axis))
    steps = np.nanmean(diff_arr, axis=axis)
    if np.isnan(steps).all():
        return 0.0

    return float(np.nanmax(steps))


def convolution(reso_params, energy_rez_factor=1 / 5, max_step=100):
    """Perform the convolution
    The maxium sampling box size in Q is (max_step, max_step ,max_step)

    Note:
        Increase the accuracy by decresing energy_rez_factor and incresing max_step
    """
    device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    # ----------------------------------------------------
    # return np.nan if repo_params is None
    # ----------------------------------------------------
    if reso_params is None:
        return np.nan
    # ----------------------------------------------------
    # calculate resolution matrix for all points
    # ----------------------------------------------------
    (qh, qk, ql), en, r0, mat = reso_params
    print(f"Calculating (Q1, Q2, Q3, E) = ({qh:.2f}, {qk:.2f}, {ql:.2f}, {en:.2f})")
    mat_hkl = quadric_proj(mat, 3)
    # ----------------------------------------------------
    # calculate the incoherent sigmas for all Q and E directions
    # ----------------------------------------------------
    sigma_qs = incoh_sigma_qs(mat_hkl)
    sigma_en_incoh = incoh_sigma_en(mat)
    num_of_sigmas = 3
    min_en, max_en = en - num_of_sigmas * sigma_en_incoh, en + num_of_sigmas * sigma_en_incoh
    sigma_en_coh = coh_sigma(mat, 3)
    # define the energy resolution to be 1/5 of the coherent sigma in energy
    en_rez = sigma_en_coh * energy_rez_factor
    # ----------------------------------------------------
    # Calculate elemental volume
    # ----------------------------------------------------
    eigenvalues = np.linalg.eigvalsh(mat_hkl)
    eval_inv_sqrt = 1 / np.sqrt(eigenvalues)
    elem_vols = np.prod(eval_inv_sqrt) * (2 * num_of_sigmas) ** 3
    # ----------------------------------------------------
    # First round, coarse grid
    # ----------------------------------------------------
    pts = [10, 10, 10]
    (vqh, vqk, vql), idx = generate_pts(sigma_qs, mat_hkl, num_of_sigmas, tuple(pts))
    disp = model_disp(vqh + qh, vqk + qk, vql + ql)
    num_bands, num_pts = disp.shape

    # Retrun zero if all dispersion is outside the relevant energy window
    if np.max(disp) < min_en or np.min(disp) > max_en:
        return 0.0
    # ----------------------------------------------------
    # determine if sampled enough based on steps along energy
    # ----------------------------------------------------
    vq = np.array((vqh, vqk, vql))  # shape: (3, num_pts)
    vqe = np.empty((4, num_bands, num_pts))
    vqe[0:3] = vq[:, None, :]
    vqe[3] = disp - en
    vqe_torch = torch.tensor(vqe, dtype=torch.float32, device=device)
    mat_torch = torch.tensor(mat, dtype=torch.float32, device=device)
    weights = compute_weights_torch(vqe_torch, mat_torch).cpu().numpy()
    # weights = compute_weights(vqe, mat)  # shape: (num_bands, num_pts)
    # Return zero if everything is outside the 5-sigma volume
    if np.min(weights) > 5**3:
        return 0.0

    # ----------------------------------------------------
    # determine Q steps based on energy steps
    # ----------------------------------------------------
    disp_arr = np.full(shape=(num_bands,) + idx.shape, fill_value=np.nan)
    disp_arr[(slice(None),) + np.nonzero(idx)] = disp
    # Compute max energy steps
    steps = [get_max_step(disp_arr, axis=i) for i in (1, 2, 3)]
    # limit the maximum in case the dispersion is too steep
    for i, (step, pt) in enumerate(zip(steps, pts)):
        if step > en_rez:
            factor = step / en_rez
            pts[i] = int(np.min((pt * factor, max_step)))

    # ----------------------------------------------------
    # Enough sampled. Calculate weight from resolution function
    # ----------------------------------------------------
    (vqh, vqk, vql), idx = generate_pts(sigma_qs, mat_hkl, num_of_sigmas, tuple(pts))
    disp = model_disp(vqh + qh, vqk + qk, vql + ql)
    _, num_pts = disp.shape

    vq = np.array((vqh, vqk, vql))  # shape: (3, num_pts)
    vqe = np.empty((4, num_bands, num_pts))
    vqe[0:3] = vq[:, None, :]
    vqe[3] = disp - en

    vqe_torch = torch.tensor(vqe, dtype=torch.float32, device=device)
    mat_torch = torch.tensor(mat, dtype=torch.float32, device=device)
    weights = compute_weights_torch(vqe_torch, mat_torch).cpu().numpy()
    # weights = compute_weights(vqe, mat)  # shape: (num_bands, num_pts)
    # ----------------------------------------------------
    # Keep only the points within the 4D ellipsoid
    # ----------------------------------------------------
    idx_keep = np.any(weights < 5**3, axis=0)
    vq_filtered = vq[:, idx_keep]
    num_pts_keep = np.count_nonzero(idx_keep)
    percent_kep = num_pts_keep / np.prod(pts) * 100
    print(f"Number of pts inside the ellipsoid = {num_pts_keep}, percentage ={percent_kep:.3f}%")

    weights_filtered = np.exp(-weights[:, idx_keep] / 2)
    inten = model_inten(*vq_filtered)
    # normalization by elementary volume size
    elem_vols /= np.prod(pts)
    det = np.linalg.det(mat)
    inten_sum = np.sum(inten * weights_filtered) * elem_vols
    return r0 * inten_sum * np.sqrt(det) / (2 * np.pi) ** 2


if __name__ == "__main__":
    device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    print(f"Using device: {device}")
    # ----------------------------------------------------
    # points being measured
    # qe_mesh has the dimension (4, n_pts_of_measurement)
    # flatten for meshed measurement
    # ----------------------------------------------------
    q1_min, q1_max, q1_step = 2, 3, 0.02
    en_min, en_max, en_step = -3, 25, 0.5
    q2 = 0
    q3 = 0

    q1 = np.linspace(q1_min, q1_max, int((q1_max - q1_min) / q1_step) + 1)
    en = np.linspace(en_min, en_max, int((en_max - en_min) / en_step) + 1)

    # calculate resolution
    vq1, vq2, vq3 = np.meshgrid(q1, q2, q3, indexing="ij")
    q_list = np.stack((vq1.ravel(), vq2.ravel(), vq3.ravel()), axis=-1)
    reso_params = resolution_matrix(hkl=q_list, en=en)

    t0 = time()
    # ------------------- multiprocessing ------------------
    # num_worker = 8
    # with ProcessPoolExecutor(max_workers=num_worker) as executor:
    #     results = executor.map(convolution, reso_params)
    # measurement_inten = np.asarray(list(results))
    # ------------------- single core ------------------
    sz = len(reso_params)
    measurement_inten = np.empty(shape=sz)
    for i in range(sz):
        measurement_inten[i] = convolution(reso_params[i])
    # --------------------------------------------------

    print(f"Convolution completed in {(t1 := time()) - t0:.4f} s")
    # total intensity should be close to S/2 *(q1_max - q1_min) * 2p*i
    total_intent = np.sum(measurement_inten) * q1_step * en_step / (q1_max - q1_min)

    # ----------------------------------------------------
    # plot 2D contour
    # ----------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 6))
    vq1, ven = np.meshgrid(q1, en, indexing="ij")
    img = ax.pcolormesh(vq1, ven, measurement_inten.reshape(np.shape(vq1)), cmap="turbo", vmin=0, vmax=0.5)

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
        + f"\n3D Convolution for {len(q1) * len(en)} points, "
        + f"completed in {t1 - t0:.3f} s"
        # + " with {num_worker:1d} cores"
    )

    plt.show()
