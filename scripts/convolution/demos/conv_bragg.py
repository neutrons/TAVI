import math
from time import time

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse


def gaussian(x, mu, sigma):
    "1D normalized Gaussian"
    return np.exp(-((x - mu) ** 2) / 2 / sigma**2) / (sigma * np.sqrt(2 * np.pi))


def bragg_peak(q1, q2, q3, en, cen=(0, 0, 0, 0), sigma=(0.05, 0.05, 0.05, 0.1), amp=1):
    """Bragg peak with cen and sigma, amp is intensity"""
    cen_q1, cen_q2, cen_q3, cen_en = cen
    sigma_q1, sigma_q2, sigma_q3, sigma_en = sigma
    return (
        amp
        * np.exp(-((q1 - cen_q1) ** 2) / 2 / sigma_q1**2)
        * np.exp(-((q2 - cen_q2) ** 2) / 2 / sigma_q2**2)
        * np.exp(-((q3 - cen_q3) ** 2) / 2 / sigma_q3**2)
        * np.exp(-((en - cen_en) ** 2) / 2 / sigma_en**2)
        / sigma_q1
        / sigma_q2
        / sigma_q3
        / sigma_en
        / (2 * np.pi) ** 2
    )


def model(vq1, vq2, vq3, ven):
    "Two splitted Bragg peaks"
    return bragg_peak(vq1, vq2, vq3, ven, cen=(-0.5, 0, 0, 0)) + bragg_peak(vq1, vq2, vq3, ven, cen=(0.5, 0, 0, 0))


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

    sigma1, sigma2 = 1, 0.2
    angle = -70
    mat = np.array(
        [
            [1 / sigma1**2, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1 / sigma2**2],
        ]
    )
    rez_mat = rotation_matrix_4d(angle).T @ mat @ rotation_matrix_4d(angle)
    r0 = 1
    return np.broadcast_to(r0, sz), np.broadcast_to(rez_mat, sz + (4, 4))


def resolution_fnc(q, m, r0):
    dim = np.size(q)
    return r0 * np.sqrt(np.linalg.det(m)) * np.exp(-q.T @ m @ q / 2) / (2 * np.pi) ** (dim / 2)


def rez_conv_4d(vq1, vq2, vq3, ven, n_sample=100000):
    t0 = time()

    r, rez_mat = resolution_matrix(vq1, vq2, vq3, ven)
    cov = np.linalg.inv(rez_mat)
    print(f"Done inverseing matriices {(t1:=time())-t0:.4f} s")

    pts_norm = np.random.default_rng().multivariate_normal(mean=(0, 0, 0, 0), cov=np.diag((1, 1, 1, 1)), size=n_sample)
    print(f"Done sampleing in {(t2:=time())-t1:.4f} s")

    eigenvalues, eigenvectors = np.linalg.eig(cov)
    pts = (
        pts_norm
        @ eigenvectors
        @ (np.sqrt(eigenvalues)[:, :, :, :, np.newaxis] * np.diag((1, 1, 1, 1)))
        @ np.transpose(eigenvectors, (0, 1, 2, 3, 5, 4))
    )
    print(f"Done transformation in {(t3:=time())-t2:.4f} s")

    pts[:, :, :, :, :, 0] += vq1[:, :, :, :, np.newaxis]
    pts[:, :, :, :, :, 1] += vq2[:, :, :, :, np.newaxis]
    pts[:, :, :, :, :, 2] += vq3[:, :, :, :, np.newaxis]
    pts[:, :, :, :, :, 3] += ven[:, :, :, :, np.newaxis]
    print(f"Done recentering in {(t4:=time())-t3:.4f} s")

    # ---------------------- direct sampling --------------------
    # sz = vq1.shape
    # pts = []
    # it = np.nditer(vq1, flags=["multi_index"])
    # for q1 in it:
    #     idx = it.multi_index
    #     pts.append(
    #         np.random.default_rng().multivariate_normal(
    #             mean=(q1, vq2[idx], vq3[idx], ven[idx]), cov=cov[idx], size=n_sample
    #         )
    #     )
    # print(f"Done recentering in {(t4:=time()-t1):.4f} s")
    # ----------------------------------------------------------

    model_cal = model(pts[:, :, :, :, :, 0], pts[:, :, :, :, :, 1], pts[:, :, :, :, :, 2], pts[:, :, :, :, :, 3])
    print(f"Done calculating intensity of the model in {(t5:=time())-t4:.4f} s")

    rez_conv = np.mean(model_cal, axis=-1)
    print(f"Done averaging in {(t6:=time())-t5:.4f} s")
    return rez_conv


if __name__ == "__main__":
    xmin, xmax = -1, 1
    ymin, ymax = -5, 5
    # ------------------- plot signal ------------------------
    nx, ny = 1000, 500
    q1 = np.linspace(xmin, xmax, nx + 1)
    en = np.linspace(ymin, ymax, ny + 1)
    q2 = 0
    q3 = 0
    # vq1, vq2, vq3, ve = np.meshgrid(q1, q2, q3, en, indexing="ij")

    fig, axes = plt.subplots(ncols=2, nrows=2, sharex=True)

    ax0 = axes[0, 0]
    vq1, ve = np.meshgrid(q1, en, indexing="ij")

    model_inten = model(*np.meshgrid(q1, q2, q3, en, indexing="ij"))
    integ_inten = np.sum(model_inten) * (xmax - xmin) / nx * (ymax - ymin) / ny
    print(f"integrated intensity of model = {integ_inten}")

    img0 = ax0.pcolormesh(vq1, ve, np.squeeze(model_inten), cmap="turbo", vmin=0)
    ax0.set_title("model")
    ax0.grid(alpha=0.6)
    # ax0.set_xlabel("Q")
    ax0.set_ylabel("E")
    fig.colorbar(img0, ax=ax0)

    ax0.add_artist(
        Ellipse(
            xy=(-0.5, 0),
            width=1 * 2,
            height=0.2 * 2,
            angle=70,
            edgecolor="w",
            facecolor="none",
            label="1-sigma ellipse",
        )
    )
    # ax0.axis("equal")
    # ------------------- plot signal cut------------------------
    ax2 = axes[0, 1]
    ax2.plot(q1, np.squeeze(model(*np.meshgrid(q1, 0, 0, 0, indexing="ij"))), "-o")
    # ax2.set_xlabel("Q")
    ax2.set_ylabel("Intensity")
    ax2.grid(alpha=0.6)
    ax2.set_title("model E = 0")

    # ------------------- plot measurement ------------------------
    xstep, ystep = 0.05, 0.5
    nx, ny = math.ceil((xmax - xmin) / xstep), math.ceil((ymax - ymin) / ystep)
    q1 = np.linspace(xmin, xmax, nx + 1)
    en = np.linspace(ymin, ymax, ny + 1)
    # vq1, vq2, vq3, ve = np.meshgrid(q1, q2, q3, en, indexing="ij")
    vq1, ve = np.meshgrid(q1, en, indexing="ij")

    ax1 = axes[1, 0]
    t_start = time()
    n_sample = 300000
    measurement_inten = rez_conv_4d(*np.meshgrid(q1, q2, q3, en, indexing="ij"), n_sample)
    t_tot = time() - t_start
    print(f"Total time = {t_tot:.4f} s")
    plt.suptitle(f"Gaussian sampling, n_sample={n_sample}, time={t_tot:.2f}")

    print(f"integrated intensity of model = {integ_inten}")
    integ_inten = np.sum(measurement_inten) * xstep * ystep
    print(f"integrated intensity of measurement = {integ_inten}")
    img1 = ax1.pcolormesh(vq1, ve, np.squeeze(measurement_inten), cmap="turbo", vmin=0)
    ax1.set_title("resolution convoluted")
    ax1.grid(alpha=0.6)
    ax1.set_xlabel("Q1")
    ax1.set_ylabel("E")
    fig.colorbar(img1, ax=ax1)

    # ------------------- plot measurement cut ------------------------
    ax3 = axes[1, 1]
    ax3.plot(q1, np.squeeze(rez_conv_4d(*np.meshgrid(q1, 0, 0, 0, indexing="ij"))), "-o")
    ax3.set_xlabel("Q1")
    ax3.set_ylabel("Intensity")
    ax3.grid(alpha=0.6)
    ax3.set_title("resolution convoluted E = 0")

    plt.tight_layout()
    plt.show()
