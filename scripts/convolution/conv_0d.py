import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse


def signal_0d():
    return 1


def model():
    "Two splitted Bragg peaks in the form of (h, k, l, intensity)"
    return np.array(
        (
            (0.5, 0, 0, 1),
            (-0.5, 0, 0, 1),
        )
    )


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
                xy=(-0.5, 0),
                width=1 * 2 * (i + 1),
                height=0.2 * 2 * (i + 1),
                angle=70,
                edgecolor="w",
                facecolor="none",
                label=f"{i+1}-sigma",
            )
        )
        ax.add_artist(
            Ellipse(
                xy=(0.5, 0),
                width=1 * 2 * (i + 1),
                height=0.2 * 2 * (i + 1),
                angle=70,
                edgecolor="w",
                facecolor="none",
                # label=f"{i+1}-sigma",
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
    r0, mat = resolution_matrix(*qe_mesh)
    # r0, mat = resolution_matrix(vq1, vq2, vq3, ven)
    det = np.linalg.det(mat)
    # points to calculate from the model
    bragg_peaks = model()

    # convolution
    exp_inten = np.zeros_like(r0, dtype=float)
    for bragg_peak in bragg_peaks:
        qh, qk, ql, inten = bragg_peak
        vqe = np.stack((qe_mesh[0] - qh, qe_mesh[1] - qk, qe_mesh[2] - ql, qe_mesh[3]), axis=1)
        prod = np.squeeze(vqe[:, np.newaxis, :] @ mat @ vqe[:, :, np.newaxis])
        exp_inten += np.exp(-prod) * inten
    exp_inten *= det * r0 / (2 * np.pi) ** 2
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
