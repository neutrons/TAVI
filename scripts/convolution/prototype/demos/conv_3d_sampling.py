import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse


def model(vq1, vq2, vq3):
    """return (energy, intensity) for given Q points
    3d FM J=-1 meV S=1"""
    sj = 1
    gamma_q = (np.cos(2 * np.pi * vq1) + np.cos(2 * np.pi * vq2) + np.cos(2 * np.pi * vq3)) / 3
    return (6 * sj * (1 - gamma_q), np.ones_like(vq1, dtype=float) / 2)


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


def plot_rez_ellipses(ax):
    for i in range(3):
        ax.add_artist(
            Ellipse(
                xy=(0, 0),
                width=1 * 2 * (i + 1),
                height=0.2 * 2 * (i + 1),
                angle=70,
                edgecolor="w",
                facecolor="none",
                label=f"{i + 1}-sigma",
            )
        )


if __name__ == "__main__":
    # points being measured
    q1_min, q1_max, q1_step = -1, 1, 0.02
    en_min, en_max, en_step = -5, 5, 0.2
    q2 = 0
    q3 = 0

    q1 = np.linspace(q1_min, q1_max, int((q1_max - q1_min) / q1_step) + 1)
    en = np.linspace(en_min, en_max, int((en_max - en_min) / en_step) + 1)
    vq1, vq2, vq3, ven = np.meshgrid(q1, q2, q3, en, indexing="ij")

    # keep the dimension then flatten
    sz = np.shape(vq1)
    qe_mesh = np.array([np.ravel(v) for v in (vq1, vq2, vq3, ven)])

    # calculate resolution matrix for all points
    # r0 has dimension of (n_q1, n_q2, n_q3, n_en)
    # mat has dimension of (n_q1, n_q2, n_q3, n_en, 4, 4)
    r0, rez_mat = resolution_matrix(*qe_mesh)
    cov = np.linalg.inv(rez_mat)

    # generate points in a sphere
    n_sample = 270_000
    n = int(n_sample ** (1 / 4))
    x = np.linspace(-3, 3, n)
    v1, v2, v3, v4 = np.meshgrid(x, x, x, x, indexing="ij")
    g = np.exp(-(v1**2) / 2) * np.exp(-(v2**2) / 2) * np.exp(-(v3**2) / 2) * np.exp(-(v4**2) / 2) / (2 * np.pi) ** 2
    sigma_cutoff = 2.5  # cut the corners beyond 2.5*sigma
    idx = g > np.max(g) * np.exp(-(sigma_cutoff**2) / 2)
    pts_norm = np.stack((v1[idx], v2[idx], v3[idx], v4[idx]), axis=1)

    weight = g[idx]
    weight /= np.sum(weight)

    # similarity transformation into ellipsoids
    eigenvalues, eigenvectors = np.linalg.eig(cov)
    mat = (
        eigenvectors
        @ (np.sqrt(eigenvalues)[:, np.newaxis] * np.diag((1, 1, 1, 1)))
        @ np.swapaxes(eigenvectors, 1, 2)  # np.transpose(eigenvectors, (0, 2, 1))
    )
    # pts has dimension (n_measurement, n_sampling, 4)
    pts = pts_norm @ mat

    fm_disp, fm_intent = model(pts[:, :, 0], pts[:, :, 1], pts[:, :, 2])
    pass

    # convolution
    exp_inten = np.zeros_like(r0, dtype=float)

    vqe = np.stack((qe_mesh[0] - qh, qe_mesh[1] - qk, qe_mesh[2] - ql, qe_mesh[3]), axis=1)
    prod = np.squeeze(vqe[:, np.newaxis, :] @ mat @ vqe[:, :, np.newaxis])
    exp_inten += np.exp(-prod / 2) * inten
    det = np.linalg.det(mat)
    exp_inten *= np.sqrt(det) * r0 / (2 * np.pi) ** 2

    exp_inten = exp_inten.reshape(sz)

    # plot 2D contour
    fig, ax = plt.subplots()
    idx = np.s_[:, 0, 0, :]
    img = ax.pcolormesh(vq1[idx], ven[idx], exp_inten[idx], cmap="turbo", vmin=0, vmax=1)
    # ax.set_title("model")
    ax.grid(alpha=0.6)
    ax.set_xlabel("Q1")
    ax.set_ylabel("E")
    ax.set_xlim((-1, 1))
    ax.set_ylim((-5, 5))
    fig.colorbar(img, ax=ax)
    plot_rez_ellipses(ax)
    ax.legend()

    plt.tight_layout()
    plt.show()
