# from concurrent.futures import ProcessPoolExecutor
import functools
from concurrent.futures import ProcessPoolExecutor
from functools import partial

import numpy as np
from numba import njit, prange


def quadric_proj(quadric: np.ndarray, idx: int) -> np.ndarray:
    """projects along one axis of the quadric

    # comparing with simple projection
    # rank = len(quadric)
    # vec /= np.sqrt(np.dot(vec, vec))
    # proj_op = np.outer(vec, vec)
    # ortho_proj = np.dot((np.identity(rank) - proj_op), quadric)"""

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


def coh_sigma(mat, axis):
    """Coherent sigma"""
    idx = int(axis)

    return 1 / np.sqrt(np.abs(mat[idx, idx]))


@njit(parallel=True, nogil=True)
def compute_weights(vqe: np.ndarray, mat: np.ndarray) -> np.ndarray:
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


@functools.cache
def generate_meshgrid(num_of_sigmas=3, num_pts=(10, 10, 10)):
    pts_qh, pts_qk, pts_ql = num_pts
    qh = np.linspace(-num_of_sigmas, num_of_sigmas, pts_qh + 1)
    qk = np.linspace(-num_of_sigmas, num_of_sigmas, pts_qk + 1)
    ql = np.linspace(-num_of_sigmas, num_of_sigmas, pts_ql + 1)
    return np.meshgrid(qh, qk, ql, indexing="ij")  # shape (3, N1, N2, N3)


def generate_pts(sigma_qs, mat_hkl, num_of_sigmas=3, num_pts=(10, 10, 10)):
    (sigma_qh_incoh, sigma_qk_incoh, sigma_ql_incoh) = sigma_qs

    vq_h, vq_k, vq_l = generate_meshgrid(num_of_sigmas, num_pts)
    vq = (vq_h * sigma_qh_incoh, vq_k * sigma_qk_incoh, vq_l * sigma_ql_incoh)

    # -------- cut the corners based on distance --------
    r_sq = np.einsum("i...,ij,j...->...", vq, mat_hkl, vq)
    idx = r_sq < num_of_sigmas**2  # Ellipsoid mask
    return (vq[0][idx], vq[1][idx], vq[2][idx]), idx
    # return vq, idx


def get_max_step(arr, axis: int):
    """Get max step along a given axis. Return zero if all NaN"""
    # shape of arr (num_bands, N1, N2, N3)
    diff_arr = np.abs(np.diff(arr, axis=axis))
    steps = np.nanmean(diff_arr, axis=axis)
    if np.isnan(steps).all():
        return 0.0

    return float(np.nanmax(steps))


def convolution(reso_params, model_disp, model_inten, energy_rez_factor=1 / 5, max_step=100):
    if reso_params is None:
        return np.nan
    # ----------------------------------------------------
    # calculate resolution matrix for all points
    # ----------------------------------------------------
    (qh, qk, ql), en, r0, mat = reso_params
    mat_hkl = quadric_proj(mat, 3)
    # ----------------------------------------------------
    # calculate the incoherent sigmas for all Q and E directions
    # ----------------------------------------------------
    # sigma_qs = [incoh_sigma_q(mat_hkl, i) for i in range(3)]
    sigma_qs = incoh_sigma_qs(mat_hkl)
    sigma_en_incoh = incoh_sigma_en(mat)
    num_of_sigmas = 3
    min_en, max_en = en - num_of_sigmas * sigma_en_incoh, en + num_of_sigmas * sigma_en_incoh
    sigma_en_coh = coh_sigma(mat, 3)
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
    # vqh, vqk, vql = vqh[idx], vqk[idx], vql[idx]
    # ----------------------------------------------------
    # determine if sampled enough based on steps along energy
    # ----------------------------------------------------
    disp = model_disp(vqh + qh, vqk + qk, vql + ql)
    num_bands, num_pts = disp.shape

    # Skip if all dispersion is outside the relevant energy window
    if np.max(disp) < min_en or np.min(disp) > max_en:
        return 0.0

    vq = np.array((vqh, vqk, vql))  # shape: (3, num_pts)
    vqe = np.empty((4, num_bands, num_pts))
    vqe[0:3] = vq[:, None, :]
    vqe[3] = disp - en

    weights = compute_weights(vqe, mat)  # shape: (num_bands, num_pts)
    # Skip if everything is outside the 5-sigma volume
    if np.min(weights) > 5**3:
        return 0.0

    # ----------------------------------------------------
    # determine Q steps based on energy steps
    # ----------------------------------------------------
    disp_arr = np.full(shape=(num_bands,) + idx.shape, fill_value=np.nan)
    disp_arr[(slice(None),) + np.nonzero(idx)] = disp

    # Compute max energy steps
    steps = [get_max_step(disp_arr, axis=i) for i in (1, 2, 3)]

    # max_step = 100  # limit the maximum in case the dispersion is too steep
    for i, (step, pt) in enumerate(zip(steps, pts)):
        if step > en_rez:
            factor = step / en_rez
            pts[i] = int(np.min((pt * factor, max_step)))

    # ----------------------------------------------------
    # Enough sampled. Calculate weight from resolution function
    # ----------------------------------------------------
    (vqh, vqk, vql), idx = generate_pts(sigma_qs, mat_hkl, num_of_sigmas, tuple(pts))
    # vqh, vqk, vql = vqh[idx], vqk[idx], vql[idx]
    disp = model_disp(vqh + qh, vqk + qk, vql + ql)
    _, num_pts = disp.shape

    # ----------------------------------------------------
    # determine Q steps based on energy steps
    # ----------------------------------------------------
    disp_arr = np.full(shape=(num_bands,) + idx.shape, fill_value=np.nan)
    disp_arr[(slice(None),) + np.nonzero(idx)] = disp

    # Compute max energy steps
    # steps = [get_max_step(disp_arr, axis=i) for i in (1, 2, 3)]
    print(f"Calculating (Q1, Q2, Q3, E) = ({qh:.2f}, {qk:.2f}, {ql:.2f}, {en:.2f})")
    # print(f"number of points = {pts}")
    # print(f"steps in energy = ({steps[0]:.2f}, {steps[1]:.2f}, {steps[2]:.2f})")

    vq = np.array((vqh, vqk, vql))  # shape: (3, num_pts)
    vqe = np.empty((4, num_bands, num_pts))
    vqe[0:3] = vq[:, None, :]
    vqe[3] = disp - en

    weights = compute_weights(vqe, mat)  # shape: (num_bands, num_pts)

    idx_keep = np.any(weights < 5**3, axis=0)
    vq_filtered = vq[:, idx_keep]
    # num_pts_keep = np.count_nonzero(idx_keep)
    # percent_kep = num_pts_keep / np.prod(pts) * 100
    # print(f"Number of pts inside the ellipsoid = {num_pts_keep}, percentage ={percent_kep:.3f}%")
    weights_filtered = np.exp(-weights[:, idx_keep] / 2)
    inten = model_inten(*vq_filtered)

    # normalization
    elem_vols /= np.prod(pts)
    det = np.linalg.det(mat)
    inten_sum = np.sum(inten * weights_filtered) * elem_vols
    return r0 * inten_sum * np.sqrt(det) / (2 * np.pi) ** 2


def calculate(model_disp, model_inten, reso_params, num_worker=8):
    conv_model = partial(convolution, model_disp=model_disp, model_inten=model_inten)

    with ProcessPoolExecutor(max_workers=num_worker) as executor:
        results = list(executor.map(conv_model, reso_params))
    return np.asarray(results)
