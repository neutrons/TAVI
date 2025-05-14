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

    for i in (1, 0):
        if not i == idx:
            mat = quadric_proj(mat, i)

    return 1 / np.sqrt(np.abs(mat[0, 0]))


def model_disp(vq1):
    """return energy for given Q points
    1d FM J=-1 meV S=1, en=2*S*J*(1-cos(Q))
    two splitted bands
    """
    sj = 1
    gamma_q = np.cos(2 * np.pi * vq1)
    disp = 2 * sj * (1 - gamma_q)
    disp = np.array((disp - 2, disp + 2))

    # reshape if only one band
    num_disp = len(disp.shape)
    if num_disp == 1:
        disp = np.reshape(disp, (1, np.size(disp)))
    return disp


def model_inten(vq1):
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


def resolution_matrix(qx0, en0):
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
                [c, -s],
                [s, c],
            ]
        )

    sigma1, sigma2 = 0.3, 0.02
    angle = -80
    mat = np.array(
        [
            [1 / sigma1**2, 0],
            [0, 1 / sigma2**2],
        ]
    )
    rez_mat = rotation_matrix_4d(angle).T @ mat @ rotation_matrix_4d(angle)
    r0 = 1
    return np.broadcast_to(r0, sz), np.broadcast_to(rez_mat, sz + (2, 2))


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


def convolutioin(qh, en):
    # ----------------------------------------------------
    # calculate resolution matrix for all points
    # ----------------------------------------------------
    r0, mat = resolution_matrix(qh, en)
    mat_hkl = quadric_proj(mat, -1)
    # ----------------------------------------------------
    # calculate the incoherent sigmas for all Q and E directions
    # ----------------------------------------------------
    sigma_qh = incoh_sigma(mat, 0)
    sigma_en = incoh_sigma(mat, 1)

    num_of_sigmas = 3
    min_qh, max_qh = qh - num_of_sigmas * sigma_qh, qh + num_of_sigmas * sigma_qh
    min_en, max_en = en - num_of_sigmas * sigma_en, en + num_of_sigmas * sigma_en
    # ----------------------------------------------------
    # determine if the dispersion energy is in the ellipsoid in a coarse grid
    # ----------------------------------------------------
    pts_q = 10
    step_qh = (max_qh - min_qh) / pts_q
    list_qh = np.linspace(min_qh, max_qh, pts_q)
    mesh_q = list_qh
    vqh = mesh_q
    # calculate dispersion on a coarse grid
    disp = model_disp(vqh)
    # get rid of the ones that are not in the ellipsoid
    max_disp, min_disp = np.max(disp), np.min(disp)
    if max_disp < min_en or min_disp > max_en:
        return 0.0  # zero intensity
    # calculate weight from resolution function

    weights_all = []
    for i in range(disp.shape[0]):
        vqe = np.stack((vqh - qh, disp[i] - en), axis=1)
        prod = np.einsum("ij,jk,ik->i", vqe, mat, vqe)
        weights_all.append(np.exp(-prod / 2))
    weights = np.max(weights_all, axis=0)
    # don't bother if the weight is already too small
    if np.max(weights) < 1e-6:
        return 0.0  # zero intensity

    # ----------------------------------------------------
    # dispersion energy is within the ellipsoid
    # ----------------------------------------------------
    # create 1D mesh
    pts_q = 10  # start with 30 points
    sampled_enough = False
    while not sampled_enough:

        step_qh = (max_qh - min_qh) / pts_q
        vqh = np.linspace(min_qh, max_qh, pts_q)
        # ----------------------------------------------------
        # calculate weight based on the dispersion
        # ----------------------------------------------------
        disp = model_disp(vqh)  # size of (n_bands, pts_q^3)

        step_en_all = []
        for i in range(disp.shape[0]):
            step_en_all.append(np.mean(np.abs(np.diff(disp[i]))))
        step_en = np.max(step_en_all, axis=0)

        if step_en > sigma_en / 3:
            sampled_enough = False
            pts_q *= 2
        else:
            sampled_enough = True

    # ----------------------------------------------------
    # Enough sampled. Calculate weight from resolution function
    # ----------------------------------------------------
    # shape of (n_bands, pts_q^3, 2)
    vqe = np.stack((np.broadcast_to(vqh, np.shape(disp)) - qh, disp - en), axis=-1)
    prod = np.einsum("ijk,kl,ijl->ij", vqe, mat, vqe)  # shape of (n_bands, pts_q^3)
    weights = np.exp(-prod / 2)
    # ----------------------------------------------------
    # trim the corners
    # ----------------------------------------------------
    cut_off = 1e-2
    idx_all = weights > cut_off
    idx = np.bitwise_or.reduce(idx_all, axis=0)
    # percent = (np.size(idx) - np.count_nonzero(idx)) / np.size(idx) * 100
    # print(f"{percent:.2f}% of points discarded.")

    # all small weights because dispersion parallel to ellipsoid
    if not np.any(idx_all):
        return 0.0  # zero intensity

    # need a correction to enforce the normalization to one
    prod = (vqh - qh) * mat_hkl * (vqh - qh)
    g = np.exp(-prod / 2) / np.sqrt(2 * np.pi) * np.sqrt(np.linalg.det(mat_hkl))
    correction = np.sum(g) * step_qh

    vqh_filtered = vqh[idx]
    inten = model_inten(vqh_filtered)
    weights_filtered = weights[:, idx] / correction
    # normalization
    det = np.linalg.det(mat)
    inten_sum = np.sum(inten * weights_filtered) * step_qh
    return r0 * inten_sum * np.sqrt(det) / (2 * np.pi)


if __name__ == "__main__":
    # ----------------------------------------------------
    # points being measured
    # qe_mesh has the dimension (4, n_pts_of_measurement)
    # flatten for meshed measurement
    # ----------------------------------------------------
    q1_min, q1_max, q1_step = -1, 1, 0.02
    en_min, en_max, en_step = -3, 7, 0.2

    q1 = np.linspace(q1_min, q1_max, int((q1_max - q1_min) / q1_step) + 1)
    en = np.linspace(en_min, en_max, int((en_max - en_min) / en_step) + 1)
    vq1, ven = np.meshgrid(q1, en, indexing="ij")

    sz = np.shape(vq1)  # sz = (n_q1, n_q2, n_q3 , n_en)
    qe_mesh = np.array([np.ravel(v) for v in (vq1, ven)])

    t0 = time()
    # ----------------------------------------------------
    # multiprocessing using concurrent future
    # ----------------------------------------------------
    num_worker = 1
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_worker) as executor:
        results = executor.map(convolutioin, *qe_mesh)
    measurement_inten = np.array(list(results)).reshape(sz)
    # measurement_inten = convolutioin(*qe_mesh).reshape(sz)
    print(f"Convolution completed in {(t1:=time())-t0:.4f} s")

    # total intensity should be close to S/2 *(q1_max - q1_min) * 2p*i
    total_intent = np.sum(measurement_inten) * q1_step * en_step / (q1_max - q1_min)
    # ----------------------------------------------------
    # plot 2D contour
    # ----------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 4))
    idx = np.s_[:, :]
    img = ax.pcolormesh(vq1[idx], ven[idx], measurement_inten[idx], cmap="turbo", vmin=0, vmax=0.5)

    ax.grid(alpha=0.6)
    ax.set_xlabel("Q1")
    ax.set_ylabel("En")
    ax.set_xlim((q1_min, q1_max))
    ax.set_ylim((en_min, en_max))

    plot_rez_ellipses(ax)
    disp = model_disp(q1)
    for i in range(np.shape(disp)[0]):
        ax.plot(q1, disp[i], "-w")
    ax.legend()
    fig.colorbar(img, ax=ax)
    ax.set_title(
        f"1D FM S=1 J=-1, total intensity = {total_intent:.3f}"
        + f"\nConvolution for {np.shape(qe_mesh)[1]} points completed in {t1-t0:.3f} s with {num_worker:1d} cores"
    )

    plt.tight_layout()
    plt.show()
