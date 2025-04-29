from time import time

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse


def bragg_peak(vq1, vq2, vq3, ven, cen=(0, 0, 0, 0), sigma=None, amp=1):
    """Bragg peak at Q=cen, intengrated intensity=amp"""

    cen_q1, cen_q2, cen_q3, cen_en = cen
    if sigma is None:  # delta-function
        inten = np.zeros_like(vq1)
        idx = np.bitwise_and(
            np.bitwise_and(vq1 >= cen_q1, vq2 >= cen_q2),
            np.bitwise_and(vq3 >= cen_q3, ven >= cen_en),
        )
        if len(idx.shape) == 2:
            idx1 = np.argmax(idx, axis=-1)
            inten[..., idx1[idx1 != 0]] = amp
            # idx0 = np.where(np.all(~idx, axis=-1))
            # inten[..., idx0] = 0
        else:
            idx1 = np.argmax(idx, axis=-1)
            if idx1 == 0 and np.all(~idx):
                pass
            else:
                inten[..., idx1] = amp
        return inten

    else:  # Gaussians
        sigma_q1, sigma_q2, sigma_q3, sigma_en = sigma
        return (
            amp
            * np.exp(-((vq1 - cen_q1) ** 2) / 2 / sigma_q1**2)
            * np.exp(-((vq2 - cen_q2) ** 2) / 2 / sigma_q2**2)
            * np.exp(-((vq3 - cen_q3) ** 2) / 2 / sigma_q3**2)
            * np.exp(-((ven - cen_en) ** 2) / 2 / sigma_en**2)
            / (sigma_q1 * sigma_q2 * sigma_q3 * sigma_en)
            / (2 * np.pi) ** 2
        )


def model(q1, q2, q3, en):
    "Two splitted Bragg peaks"
    inten1 = bragg_peak(q1, q2, q3, en, cen=(-0.5, 0, 0, 0))
    inten2 = bragg_peak(q1, q2, q3, en, cen=(0.5, 0, 0, 0))
    # inten1 = bragg_peak(
    #     q1, q2, q3, en, cen=(-0.5, 0, 0, 0), sigma=(0.05, 0.05, 0.05, 0.2)
    # )
    # inten2 = bragg_peak(
    #     q1, q2, q3, en, cen=(0.5, 0, 0, 0), sigma=(0.05, 0.05, 0.05, 0.1)
    # )
    return inten1 + inten2


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


def rez_conv_4d(vq1, vq2, vq3, ven, n_sample=300000):
    t0 = time()

    r, rez_mat = resolution_matrix(vq1, vq2, vq3, ven)
    cov = np.linalg.inv(rez_mat)
    print(f"Done inverseing matriices {(t1:=time())-t0:.4f} s")

    n = int(n_sample ** (1 / 4))
    x = np.linspace(-3, 3, n)
    v1, v2, v3, v4 = np.meshgrid(x, x, x, x, indexing="ij")
    g = (
        np.exp(-(v1**2) / 2)
        * np.exp(-(v2**2) / 2)
        * np.exp(-(v3**2) / 2)
        * np.exp(-(v4**2) / 2)
        / (2 * np.pi) ** 2
    )
    idx = g > np.max(g) / 15  # cut the corners
    pts_norm = np.stack((v1[idx], v2[idx], v3[idx], v4[idx]), axis=1)
    print(f"Done generating gaussian in {(t2:=time())-t1:.4f} s")

    eigenvalues, eigenvectors = np.linalg.eig(cov)
    mat = (
        eigenvectors
        @ (np.sqrt(eigenvalues)[:, np.newaxis] * np.diag((1, 1, 1, 1)))
        @ np.swapaxes(eigenvectors, 1, 2)  # np.transpose(eigenvectors, (0, 2, 1))
    )
    pts = pts_norm @ mat
    print(f"Done transformation in {(t3:=time())-t2:.4f} s")

    vqe = np.stack((vq1, vq2, vq3, ven), axis=1)
    pts += vqe[:, np.newaxis, :]
    print(f"Done recentering in {(t4:=time())-t3:.4f} s")

    model_cal = model(*np.unstack(pts, axis=-1))
    print(f"Done calculating intensity of the model in {(t5:=time())-t4:.4f} s")

    rez_conv = np.mean(model_cal * np.ravel(g[idx]), axis=-1)
    print(f"Done averaging in {(t6:=time())-t5:.4f} s")
    return rez_conv


def plot_rez_ellipses(ax):
    for i in range(3):
        for xy in ((-0.5, 0), (0.5, 0)):
            ax.add_artist(
                Ellipse(
                    xy=xy,
                    width=1 * 2 * (i + 1),
                    height=0.2 * 2 * (i + 1),
                    angle=70,
                    edgecolor="w",
                    facecolor="none",
                    label=f"{i+1}-sigma",
                )
            )


if __name__ == "__main__":
    q1_min, q1_max, q1_step = -1, 1, 0.02
    en_min, en_max, en_step = -5, 5, 0.1
    q2 = 0
    q3 = 0
    # ------------------- plot signal ------------------------

    fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True)

    q1 = np.linspace(q1_min, q1_max, int((q1_max - q1_min) / q1_step) + 1)
    en = np.linspace(en_min, en_max, int((en_max - en_min) / en_step) + 1)
    vq1, ven = np.meshgrid(q1, en, indexing="ij")

    qe_mesh = (np.ravel(v) for v in np.meshgrid(q1, q2, q3, en, indexing="ij"))
    model_inten = model(*qe_mesh).reshape(vq1.shape) / q1_step / en_step

    integ_inten = np.sum(model_inten) * q1_step * en_step
    print(f"integrated intensity of model = {integ_inten}")

    img0 = ax0.pcolormesh(vq1, ven, model_inten, cmap="turbo", vmin=0)
    ax0.set_title("model")
    ax0.grid(alpha=0.6)
    # ax0.set_xlabel("Q")
    ax0.set_ylabel("E")
    fig.colorbar(img0, ax=ax0)

    plot_rez_ellipses(ax0)

    # ------------------- plot measurement ------------------------
    q1_step, en_step = 0.1, 0.2
    q1 = np.linspace(q1_min, q1_max, int((q1_max - q1_min) / q1_step))
    en = np.linspace(en_min, en_max, int((en_max - en_min) / en_step))
    vq1, ve = np.meshgrid(q1, en, indexing="ij")

    t_start = time()
    n_sample = 500000
    qe_mesh = (np.ravel(v) for v in np.meshgrid(q1, q2, q3, en, indexing="ij"))
    measurement_inten = (
        rez_conv_4d(*qe_mesh, n_sample=n_sample).reshape(vq1.shape) / q1_step / en_step
    )

    integ_inten = np.sum(measurement_inten) * q1_step * en_step
    print(f"integrated intensity of measurement = {integ_inten}")
    t_tot = time() - t_start
    print(f"Total time = {t_tot:.4f} s")

    img1 = ax1.pcolormesh(vq1, ve, measurement_inten, cmap="turbo", vmin=0)
    ax1.set_title("resolution convoluted")
    ax1.grid(alpha=0.6)
    ax1.set_xlabel("Q1")
    ax1.set_ylabel("E")
    fig.colorbar(img1, ax=ax1)

    plot_rez_ellipses(ax1)
    plt.suptitle(f"Histogramming, n_sample={n_sample}, time={t_tot:.2f}")
    plt.tight_layout()

    # ----------------------------------------------------

    fig, (ax2, ax3) = plt.subplots(nrows=2, sharex=True)

    # ------------------- plot signal cut------------------------
    model_inten = model(q1=q1, q2=0, q3=0, en=0) / q1_step
    ax2.plot(q1, model_inten, "-o")
    # ax2.set_xlabel("Q")
    ax2.set_ylabel("Intensity")
    ax2.grid(alpha=0.6)
    ax2.set_title("model E = 0")

    # # ------------------- plot measurement cut ------------------------
    qe_mesh = (np.ravel(v) for v in np.meshgrid(q1, 0, 0, 0, indexing="ij"))
    measurement_inten = rez_conv_4d(*qe_mesh, n_sample).reshape(q1.shape)
    ax3.plot(q1, measurement_inten, "-o")
    ax3.set_xlabel("Q1")
    ax3.set_ylabel("Intensity")
    ax3.grid(alpha=0.6)
    ax3.set_title("resolution convoluted E = 0")

    plt.tight_layout()
    plt.show()
