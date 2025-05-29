from time import time

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse


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
    angle = -30
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
        for xy in ((-0.5, 0), (0.5, 0)):
            ax.add_artist(
                Ellipse(
                    xy=xy,
                    width=1 * 2 * (i + 1),
                    height=0.2 * 2 * (i + 1),
                    angle=70,
                    edgecolor="w",
                    facecolor="none",
                    label=f"{i + 1}-sigma",
                )
            )


if __name__ == "__main__":
    vq1, vq2, vq3, ven = (2,), (0,), (0,), (0,)
    n_sample = 270_000
    fig, axes = plt.subplots(ncols=3, figsize=(12, 4))

    t0 = time()

    r, rez_mat = resolution_matrix(vq1, vq2, vq3, ven)
    cov = np.linalg.inv(rez_mat)
    print(f"Done inverseing matriices {(t1 := time()) - t0:.4f} s")

    n = int(n_sample ** (1 / 4))
    x = np.linspace(-3, 3, n)
    v1, v2, v3, v4 = np.meshgrid(x, x, x, x, indexing="ij")
    g = np.exp(-(v1**2) / 2) * np.exp(-(v2**2) / 2) * np.exp(-(v3**2) / 2) * np.exp(-(v4**2) / 2) / (2 * np.pi) ** 2
    sigma_cutoff = 2.5  # cut the corners beyond 2.5*sigma
    idx = g > np.max(g) * np.exp(-(sigma_cutoff**2) / 2)
    pts_norm = np.stack((v1[idx], v2[idx], v3[idx], v4[idx]), axis=1)

    weight = g[idx]
    weight /= np.sum(weight)
    print(f"Done generating gaussian in {(t2 := time()) - t1:.4f} s")

    ax = axes[0]
    ax.plot(pts_norm[:, 0], pts_norm[:, 3], ".")
    ax.set_xlim((-3, 3))
    ax.set_ylim((-3, 3))
    ax.set_aspect("equal")
    ax.set_xlabel("q1")
    ax.set_ylabel("en")
    ax.grid(alpha=0.6)

    eigenvalues, eigenvectors = np.linalg.eig(cov)
    mat = (
        eigenvectors
        @ (np.sqrt(eigenvalues)[:, np.newaxis] * np.diag((1, 1, 1, 1)))
        @ np.swapaxes(eigenvectors, 1, 2)  # np.transpose(eigenvectors, (0, 2, 1))
    )
    pts = pts_norm @ mat
    # max_dist = np.max(np.linalg.norm(pts, axis=2), axis=1, keepdims=True)
    print(f"Done transformation in {(t3 := time()) - t2:.4f} s")

    ax = axes[1]
    ax.plot(pts[0, :, 0], pts[0, :, 3], ".")
    ax.set_xlim((-3, 3))
    ax.set_ylim((-3, 3))
    ax.set_aspect("equal")
    ax.set_xlabel("q1")
    ax.set_ylabel("en")
    ax.grid(alpha=0.6)
    ax.legend()

    vqe = np.stack((vq1, vq2, vq3, ven), axis=1)
    pts += vqe[:, np.newaxis, :]
    print(f"Done recentering in {(t4 := time()) - t3:.4f} s")

    ax = axes[2]
    ax.plot(pts[0, :, 0], pts[0, :, 3], ".")
    ax.set_xlim((-1, 5))
    ax.set_ylim((-3, 3))
    ax.set_aspect("equal")
    ax.set_xlabel("q1")
    ax.set_ylabel("en")
    ax.grid(alpha=0.6)

    # putting a delta-sharpe Bragg peak at (2,0,0,en=0)
    q_bragg = (2, 0, 0, 0)
    # determin if the peak is inside the ellipsoid

    pass
    # norms = np.linalg.norm(q_bragg - pts, axis=2)
    # norms_in = norms[norms < max_dist]
    # idx0 = np.argmin(norms, axis=1, keepdims=True)

    # ax.plot(pts[0, idx0[0], 0], pts[0, idx0[0], 3], "o", label=f"Bragg {q_bragg}")
    ax.legend()
    # model_cal = model(*np.unstack(pts, axis=-1))
    # print(f"Done calculating intensity of the model in {(t5:=time())-t4:.4f} s")

    # rez_conv = np.mean(model_cal * g[idx], axis=-1)
    # print(f"Done averaging in {(t6:=time())-t5:.4f} s")

    plt.suptitle("Histogramming and simmilarity transformation")
    plt.tight_layout()

    plt.show()
