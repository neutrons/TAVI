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
            mat = quadric_proj(mat, i)

    return 1 / np.sqrt(np.abs(mat[0, 0]))


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
                label=f"{i + 1}-sigma",
            )
        )


def generate_pts_simple_cubic(num_of_sigmas=3, pts_q=5):
    x = np.linspace(-num_of_sigmas, num_of_sigmas, pts_q + 1)
    v1, v2, v3 = np.meshgrid(x, x, x, indexing="ij")
    # -------- keep all points ---------
    # vs = np.asarray([np.ravel(v) for v in (v1, v2, v3)])
    # pts_norm = np.asarray(vs)
    # -------- cut the corners based on distance --------
    r_sq = v1**2 + v2**2 + v3**2
    idx = r_sq < num_of_sigmas**2
    # ----------------------------------------
    pts_norm = np.asarray((v1[idx], v2[idx], v3[idx]))
    return pts_norm


def generate_pts(pts, n_round: int, step_q: float):
    if n_round == 0:
        return pts, 1.0

    quotient, remainder = divmod(n_round + 1, 2)
    x = step_q / (2**quotient)
    match remainder:
        case 0:  # generate face centers
            print(f"shifted by ({x:.4f},{x:.4f},0)")
            pts_shifted = np.concatenate(
                (
                    pts + np.array((x, x, 0))[:, None],
                    pts + np.array((0, x, x))[:, None],
                    pts + np.array((x, 0, x))[:, None],
                ),
                axis=1,
            )
            elem_vols_fraction = 2 / (8**quotient)
        case 1:  # generate body center and edge centers
            print(f"shifted by ({-x:.4f},{-x:.4f},{-x:.4f})")
            pts_shifted = pts + np.array((-x, -x, -x))[:, None]
            elem_vols_fraction = 1 / (8**quotient)

    return pts_shifted, elem_vols_fraction


def convolution(qh, qk, ql, en):
    # ----------------------------------------------------
    # calculate resolution matrix for all points
    # ----------------------------------------------------
    r0, mat_qe = resolution_matrix(qh, qk, ql, en)
    det = np.linalg.det(mat_qe)
    mat_q = quadric_proj(mat_qe, 3)
    # ----------------------------------------------------
    # calculate the incoherent sigmas for E directions
    # ----------------------------------------------------
    sigma_en = incoh_sigma(mat_qe, 3)
    num_of_sigmas = 3
    min_en, max_en = en - num_of_sigmas * sigma_en, en + num_of_sigmas * sigma_en
    # ----------------------------------------------------
    # similarity transformation
    # ----------------------------------------------------
    eigenvalues, eigenvectors = np.linalg.eig(mat_q)
    eval_inv_sqrt = 1 / np.sqrt(eigenvalues)
    trans_mat = eigenvectors.T @ np.diag(eval_inv_sqrt) @ eigenvectors
    # ----------------------------------------------------
    # initial parameters setup
    # ----------------------------------------------------
    pts_q = 20
    pts_norm = generate_pts_simple_cubic(num_of_sigmas, pts_q)
    step_q = 2 * num_of_sigmas / pts_q
    elem_vols_init = step_q**3 * np.prod(eval_inv_sqrt)
    # _, num_pts_q = np.shape(pts_norm)

    n_round = 0
    sampled_enough = False
    print("#" * 20)
    while not sampled_enough:
        pts_norm_shifted, elem_vols_fraction = generate_pts(pts_norm, n_round, step_q)
        vqh, vqk, vql = trans_mat @ pts_norm_shifted
        # ----------------------------------------------------
        # calculate dispersion nergy
        # ----------------------------------------------------
        disp = model_disp(vqh + qh, vqk + qk, vql + ql)
        num_bands, num_pts = disp.shape
        # ----------------------------------------------------
        # get rid of the ones that are not in the ellipsoid
        # ----------------------------------------------------
        # if n_round == 0 or n_round == 1:
        max_disp, min_disp = np.max(disp), np.min(disp)
        if max_disp < min_en or min_disp > max_en:
            print("out of bounds")
            return 0.0  # zero intensity
        # ----------------------------------------------------
        # determine if sampled enough
        # ----------------------------------------------------
        vq = np.asarray((vqh, vqk, vql))  # shape: (3, num_pts)
        vqe = np.empty((4, num_bands, num_pts))
        vqe[0:3] = vq[:, None, :]
        vqe[3] = disp - en
        prod = np.einsum("ijk,il,ljk->jk", vqe, mat_qe, vqe)
        weights = np.exp(-prod / 2)
        # don't bother if the weight is already too small
        if np.max(weights) < 1e-9:
            print("max of weight < 1e-9")
            return 0.0  # zero intensity

        # ----------------------------------------------------
        # trim the corners in (Q,E)
        # ----------------------------------------------------
        cut_off = 1e-9
        idx_all = weights > cut_off
        idx = np.any(idx_all, axis=0)
        # all small weights because dispersion parallel to ellipsoid
        if not np.any(idx_all):
            print("all weights < 1e-9")
            return 0.0  # zero intensity
        vqh_filtered, vqk_filtered, vql_filtered = vq[:, idx]
        # ----------------------------------------------------
        # calculate intensity
        # ----------------------------------------------------
        elem_vols = elem_vols_init * elem_vols_fraction

        # normalization
        inten = model_inten(vqh_filtered + qh, vqk_filtered + qk, vql_filtered + ql)
        inten_sum = np.sum(inten * weights[:, idx]) * elem_vols
        inten = r0 * inten_sum * np.sqrt(det) / (2 * np.pi) ** 2

        if n_round == 0:
            n_round = 1
            idx_list = [1]
            inten_list = [inten]
            pts_norm = pts_norm_shifted
            print("#" * 20)
            print(f"n=0 for (Q,E)=({qh:.2f},{en:.2f}), n_pts={np.shape(pts_norm_shifted)[1]}")
            print(
                f"first point = ({pts_norm_shifted[0, 0]:.3f}, {pts_norm_shifted[1, 0]:.3f}, {pts_norm_shifted[2, 0]:.3f})"
            )
            print(f"volume fraction = {elem_vols_fraction}")
            continue

        last_inten = inten_list[-1]
        print("#" * 20)
        print(f"n={n_round} for (Q,E)=({qh:.2f},{en:.2f}), n_pts={np.shape(pts_norm_shifted)[1]}")
        print(
            f"first point = ({pts_norm_shifted[0, 0]:.3f}, {pts_norm_shifted[1, 0]:.3f}, {pts_norm_shifted[2, 0]:.3f})"
        )
        print(f"volume fraction = {elem_vols_fraction}")
        print(f"last intensity = {last_inten}")
        print(f"new intensity = {inten}")
        if (n_round + 1) % 2 == 0:
            inten = inten + last_inten / 4
        else:
            inten = inten + last_inten / 2
        print(f"averaged intensity = {inten}")

        if np.abs(inten - last_inten) / last_inten < 1e-1:
            sampled_enough = True
        elif np.shape(pts_norm)[1] > 1e6:
            sampled_enough = True
        else:
            pts_norm = np.concatenate((pts_norm, pts_norm_shifted), axis=1)
            n_round += 1

        # fig, ax = plt.subplots()
        # ax.plot(pts_norm[0], pts_norm[1], "o")
        # plt.show()

        idx_list.append(n_round + 1)
        inten_list.append(inten)
    print(f"intensity list = {inten_list}")
    # print(f"final intensity = {inten}")
    return inten_list[-1]


if __name__ == "__main__":
    # ----------------------------------------------------
    # points being measured
    # qe_mesh has the dimension (4, n_pts_of_measurement)
    # flatten for meshed measurement
    # ----------------------------------------------------
    q1_min, q1_max, q1_step = 0, 0.95, 0.01
    # q1_min, q1_max, q1_step = -0.26, -0.24, 0.01
    en_min, en_max, en_step = -3, 25, 0.5
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
