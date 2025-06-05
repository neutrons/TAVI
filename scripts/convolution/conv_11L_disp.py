from concurrent.futures import ProcessPoolExecutor
from time import time

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axisartist import Axes
from numba import njit, prange

from tavi.instrument.tas import TAS
from tavi.plotter import Plot2D
from tavi.sample import Sample
from tavi.ub_algorithm import (
    UBConf,
    b_mat_from_ub_matrix,
    mantid_to_spice,
    plane_normal_from_two_peaks,
    u_mat_from_ub_matrix,
    uv_to_ub_matrix,
)


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


def get_max_step(arr, axis: int):
    """Get max step along a given axis. Reutrn zero if all NaN"""
    # shape of arr (num_bands, N1, N2, N3)
    diff_arr = np.abs(np.diff(arr, axis=axis))
    steps = np.nanmean(diff_arr, axis=axis)
    # print(f"max_step={step:.3f}")
    return np.nanmax(steps, out=np.array(0.0))


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
    vq = np.meshgrid(qh_list, qk_list, ql_list, indexing="ij")  # shape (3, N1, N2, N3)
    # -------- cut the corners based on distance --------
    r_sq = np.einsum("i...,ij,j...->...", vq, mat_hkl, vq)
    idx = r_sq < num_of_sigmas**2  # Ellipsoid mask

    # return vq, idx
    return (vq[0][idx], vq[1][idx], vq[2][idx]), idx


def convolution(reso_params):
    # ----------------------------------------------------
    # calculate resolution matrix for all points
    # ----------------------------------------------------
    if reso_params is None:
        return np.nan
    (qh, qk, ql), en, r0, mat = reso_params
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
    # First round, coarse grid
    # ----------------------------------------------------
    pts = [10, 10, 10]
    (vqh, vqk, vql), idx = generate_pts(sigma_qs, num_of_sigmas, pts, mat_hkl)
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

    max_step = 100  # limit the maximum in case the dispersion is too steep
    for i, (step, pt) in enumerate(zip(steps, pts)):
        if step > en_rez:
            factor = step / en_rez
            pts[i] = int(np.min((pt * factor, max_step)))

    # ----------------------------------------------------
    # Enough sampled. Calculate weight from resolution function
    # ----------------------------------------------------
    (vqh, vqk, vql), idx = generate_pts(sigma_qs, num_of_sigmas, pts, mat_hkl)
    # vqh, vqk, vql = vqh[idx], vqk[idx], vql[idx]
    disp = model_disp(vqh + qh, vqk + qk, vql + ql)
    _, num_pts = disp.shape

    # ----------------------------------------------------
    # determine Q steps based on energy steps
    # ----------------------------------------------------
    disp_arr = np.full(shape=(num_bands,) + idx.shape, fill_value=np.nan)
    disp_arr[(slice(None),) + np.nonzero(idx)] = disp

    # Compute max energy steps
    steps = [get_max_step(disp_arr, axis=i) for i in (1, 2, 3)]
    print(f"(Q1, Q2, Q3, E) = ({qh:.2f}, {qk:.2f}, {ql:.2f}, {en:.2f})")
    print(f"number of points = {pts}")
    print(f"steps in energy = ({steps[0]:.2f}, {steps[1]:.2f}, {steps[2]:.2f})")

    vq = np.array((vqh, vqk, vql))  # shape: (3, num_pts)
    vqe = np.empty((4, num_bands, num_pts))
    vqe[0:3] = vq[:, None, :]
    vqe[3] = disp - en

    weights = compute_weights(vqe, mat)  # shape: (num_bands, num_pts)

    idx_keep = np.any(weights < 5**3, axis=0)
    vq_filtered = vq[:, idx_keep]
    num_pts_keep = np.count_nonzero(idx_keep)
    percent_kep = num_pts_keep / np.prod(pts) * 100
    print(f"Number of pts inside the ellipsoid = {num_pts_keep}, percentage ={percent_kep:.3f}%")
    weights_filtered = np.exp(-weights[:, idx_keep] / 2)
    inten = model_inten(*vq_filtered)

    # normalization
    elem_vols /= np.prod(pts)
    det = np.linalg.det(mat)
    inten_sum = np.sum(inten * weights_filtered) * elem_vols
    return r0 * inten_sum * np.sqrt(det) / (2 * np.pi) ** 2


if __name__ == "__main__":
    # setup instrument
    instrument_config_json_path = "./src/tavi/instrument/instrument_params/hb3.json"
    hb3 = TAS(fixed_ef=14.7)
    hb3.load_instrument_params_from_json(instrument_config_json_path)
    # set up a cubic sample
    lattice_params = (10, 10, 10, 90, 90, 90)
    sample = Sample(lattice_params)
    sample.set_mosaic(30, 30)
    # set up sample orientation, u along ki, v in plane, u cross v is up
    # when all goniometer angles are zeros
    u = (1, 1, 0)
    v = (0, 0, 1)
    # set up the scattering plane, (100) and (010) in two peaks in plane
    peaks = ((1, 1, 0), (0, 0, 1))
    ub_matrix_mantid = uv_to_ub_matrix(u, v, lattice_params)
    plane_normal_mantid, in_plane_ref_mantid = plane_normal_from_two_peaks(
        u_mat_from_ub_matrix(ub_matrix_mantid),
        b_mat_from_ub_matrix(ub_matrix_mantid),
        *peaks,
    )
    sample.ub_conf = UBConf(
        ub_mat=mantid_to_spice(ub_matrix_mantid),
        plane_normal=mantid_to_spice(plane_normal_mantid),
        in_plane_ref=mantid_to_spice(in_plane_ref_mantid),
    )
    hb3.mount_sample(sample)

    # ----------------------------------------------------
    # points being measured
    # ----------------------------------------------------
    q1_min, q1_max, q1_step = -3, 3, 0.02
    en_min, en_max, en_step = -3, 25, 0.5

    q1 = np.linspace(q1_min, q1_max, int((q1_max - q1_min) / q1_step) + 1)
    en = np.linspace(en_min, en_max, int((en_max - en_min) / en_step) + 1)
    q_list = np.array([(1, 1, l) for l in q1])

    reso_params = [
        (reso.hkl, reso.en, reso.r0, reso.mat) if reso is not None else None
        for reso in hb3.cooper_nathans(hkl=q_list, en=en)
    ]

    t0 = time()
    num_worker = 8
    with ProcessPoolExecutor(max_workers=num_worker) as executor:
        results = executor.map(convolution, reso_params)
    measurement_inten = np.asarray(list(results))

    print(f"Convolution completed in {(t1 := time()) - t0:.4f} s")

    # ----------------------------------------------------
    # plot 2D contour
    # ----------------------------------------------------
    # calculate and plot resolution
    q1_list = np.linspace(q1_min, q1_max, int((q1_max - q1_min) / (q1_step * 10)) + 1)
    en_list = np.linspace(en_min, en_max, int((en_max - en_min) / (en_step * 10)) + 1)
    q_list = np.array([(1, 1, l) for l in q1_list])
    rez_list = hb3.cooper_nathans(hkl=q_list, en=en_list, projection=((1, 1, 0), (-1, 1, 0), (0, 0, 1)))

    p = Plot2D()
    for rez in rez_list:
        if rez is None:
            continue
        e_co = rez.get_ellipse(axes=(2, 3), PROJECTION=False)
        e_inco = rez.get_ellipse(axes=(2, 3), PROJECTION=True)
        p.add_reso(e_co, c="w", linestyle="solid")
        p.add_reso(e_inco, c="w", linestyle="dashed")
    # create plot
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, axes_class=Axes)
    p.plot(ax)

    # overplot contour
    vq1, ven = np.meshgrid(q1, en)
    im = ax.pcolormesh(
        vq1,
        ven,
        measurement_inten.reshape(np.shape(vq1)),
        cmap="turbo",
        vmin=0,
        vmax=0.005,
    )
    fig.colorbar(im, ax=ax)

    # plot dispersion
    disp = model_disp(np.ones_like(q1), np.ones_like(q1), q1)
    for i in range(np.shape(disp)[0]):
        ax.plot(q1, disp[i], "-w")

    ax.set_xlim((q1_min, q1_max))
    ax.set_ylim((en_min, en_max))

    ax.set_title(
        "3D FM S=1 J=-5"
        + f"\n3D Convolution for {len(q1) * len(en)} points, "
        + f"completed in {t1 - t0:.3f} s with {num_worker:1d} cores"
    )
    ax.grid(alpha=0.6)
    # plt.tight_layout()
    plt.show()
