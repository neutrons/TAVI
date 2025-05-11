import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse


def signal_0d():
    return 1


def model(q1, q2, q3, en):
    "Two splitted Bragg peaks"
    pass


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
    qe_mesh = (np.ravel(v) for v in (vq1, vq2, vq3, ven))

    # calculate resolution matrix for all points
    # r0 has dimension of (n_q1 * n_q2 * n_q3 *n_en, )
    # mat has dimension of (n_q1 * n_q2 * n_q3 * n_en, 4, 4)
    r0, mat = resolution_matrix(*qe_mesh)

    # points to calculate from the model
    mq1_min, mq1_max, mq1_step = -2, 2, 1e-3
    mq2_min, mq2_max, mq2_step = -0.5, 0.5, 1e-3
    mq3_min, mq3_max, mq3_step = -0.5, 0.5, 1e-3
    men_min, men_max, men_step = -8, 8, 1e-2

    mq1 = np.linspace(mq1_min, mq1_max, int((mq1_max - mq1_min) / mq1_step) + 1)
    mq2 = np.linspace(mq2_min, mq2_max, int((mq2_max - mq2_min) / mq2_step) + 1)
    mq3 = np.linspace(mq3_min, mq3_max, int((mq3_max - mq3_min) / mq3_step) + 1)
    men = np.linspace(men_min, men_max, int((men_max - men_min) / en_step) + 1)
    vmq1, vmq2, vmq3, vmen = np.meshgrid(mq1, mq2, mq3, men, indexing="ij")

    plt.tight_layout()
    plt.show()
