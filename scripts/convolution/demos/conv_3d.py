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

    sigma1, sigma2 = 0.3, 0.02
    angle = -80
    mat = np.array(
        [
            [1 / sigma1**2, 0, 0, 0],
            [0, 1 / 0.1**2, 0, 0],
            [0, 0, 1 / 0.1**2, 0],
            [0, 0, 0, 1 / sigma2**2],
        ]
    )
    rez_mat = rotation_matrix_4d(angle).T @ mat @ rotation_matrix_4d(angle)
    r0 = 1
    return np.broadcast_to(r0, sz), np.broadcast_to(rez_mat, sz + (4, 4))


def plot_rez_ellipses(ax):
    for i in range(3):
        sigma1, sigma2 = 0.3, 0.02

        ax.add_artist(
            Ellipse(
                xy=(0, 0),
                width=sigma1 * 2 * (i + 1),
                height=sigma2 * 2 * (i + 1),
                angle=80,
                edgecolor="w",
                facecolor="none",
                label=f"{i+1}-sigma",
            )
        )


if __name__ == "__main__":
    # points being measured
    q1_min, q1_max, q1_step = -2, 2, 0.1
    en_min, en_max, en_step = -5, 12, 0.2
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

    # points to calculate from the model
    mq1_min, mq1_max, mq1_step = -2, 2, 0.05
    mq2_min, mq2_max, mq2_step = -0.5, 0.5, 0.05
    mq3_min, mq3_max, mq3_step = -0.5, 0.5, 0.05

    mq1 = np.linspace(mq1_min, mq1_max, int((mq1_max - mq1_min) / mq1_step) + 1)
    mq2 = np.linspace(mq2_min, mq2_max, int((mq2_max - mq2_min) / mq2_step) + 1)
    mq3 = np.linspace(mq3_min, mq3_max, int((mq3_max - mq3_min) / mq3_step) + 1)
    vmq1, vmq2, vmq3 = np.meshgrid(mq1, mq2, mq3, indexing="ij")

    msz = np.shape(vmq1)
    mqe_mesh = np.array([np.ravel(v) for v in (vmq1, vmq2, vmq3)])

    fm_disp, fm_intent = model(*mqe_mesh)

    # convolution
    exp_inten = np.zeros_like(r0, dtype=float)

    # for i in range(np.shape(mqe_mesh)[-1]):
    #     qh, qk, ql = mqe_mesh[:, i]
    #     vqe = np.stack((qe_mesh[0] - qh, qe_mesh[1] - qk, qe_mesh[2] - ql, qe_mesh[3] - fm_disp[i]), axis=1)
    #     prod = np.squeeze(vqe[:, np.newaxis, :] @ mat @ vqe[:, :, np.newaxis])
    #     exp_inten += np.exp(-prod / 2) * fm_intent[i]

    # Assuming mqe_mesh has shape (3, N)
    qh, qk, ql = mqe_mesh
    vqe = np.stack(
        (
            qe_mesh[0][:, None] - qh,
            qe_mesh[1][:, None] - qk,
            qe_mesh[2][:, None] - ql,
            qe_mesh[3][:, None] - fm_disp,
        ),
        axis=2,
    )  # shape: (M, N, 4)

    # Compute vqe @ mat @ vqe^T for each M and N
    # einsum does: vqe[m,n,:] @ mat[m,:,:] @ vqe[m,n,:].T -> shape: (M, N)
    prod = np.einsum("mni,mij,mnj->mn", vqe, mat, vqe)
    # Compute exponential terms and weight by fm_intent
    weighted_exp = np.exp(-prod / 2) * fm_intent  # shape: (M, N)

    # Sum across N (the original loop's dimension)
    exp_inten += np.sum(weighted_exp, axis=1)
    det = np.linalg.det(mat)
    exp_inten *= np.sqrt(det) * r0 / (2 * np.pi) ** 2

    exp_inten = exp_inten.reshape(sz)

    # plot 2D contour
    fig, ax = plt.subplots()
    idx = np.s_[:, 0, 0, :]
    img = ax.pcolormesh(vq1[idx], ven[idx], exp_inten[idx], cmap="turbo", vmin=0)
    # ax.set_title("model")
    ax.grid(alpha=0.6)
    ax.set_xlabel("Q1")
    ax.set_ylabel("E")
    ax.set_xlim((q1_min, q1_max))
    ax.set_ylim((en_min, en_max))
    fig.colorbar(img, ax=ax)
    plot_rez_ellipses(ax)
    ax.legend()

    plt.tight_layout()
    plt.show()
