import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse


def gaussian(x, mu, sigma):
    "1D normalized Gaussian"
    return np.exp(-((x - mu) ** 2) / 2 / sigma**2) / (sigma * np.sqrt(2 * np.pi))


def bragg_peak(q1, q2, q3, en, cen=(0, 0, 0, 0), sigma=(0.05, 0.05, 0.05, 0.1)):
    """Bragg peak with cen and sigma, amp is intensity"""
    cen_qx, cen_qy, cen_qz, cen_en = cen
    sigma_qx, sigma_qy, sigma_qz, sigma_en = sigma
    amp = 1

    q1, q2, q3, en = np.meshgrid(q1, q2, q3, en, indexing="ij")
    return np.squeeze(
        amp
        * gaussian(q1, cen_qx, sigma_qx)
        * gaussian(q2, cen_qy, sigma_qy)
        * gaussian(q3, cen_qz, sigma_qz)
        * gaussian(en, cen_en, sigma_en)
    )


def signal(q1, q2, q3, en):
    "Two splitted Bragg peaks"
    return bragg_peak(q1, q2, q3, en, cen=(-0.5, 0, 0, 0)) + bragg_peak(q1, q2, q3, en, cen=(0.5, 0, 0, 0))


def resolution_matrix():
    """Fake resoltuion matrix mat and prefactor r0
    r0 is a constant, rez_mat is a symmatric positive 4 by 4 matrix
    """

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
    return r0, rez_mat


def resolution_fnc(q, m, r0):
    dim = np.size(q)
    return r0 * np.sqrt(np.linalg.det(m)) * np.exp(-q.T @ m @ q / 2) / (2 * np.pi) ** (dim / 2)


def rez_conv_4d(q1, q2, q3, en):

    q1, q2, q3, en = np.meshgrid(q1, q2, q3, en, indexing="ij")
    sz = q1.shape
    result = []
    for qx0, qy0, qz0, en0 in zip(q1.flat, q2.flat, q3.flat, en.flat):
        r0, rez_mat = resolution_matrix()
        cov = np.linalg.inv(rez_mat)

        sampling_pts = 5000
        rez_conv = 0
        rng = np.random.default_rng()
        pts = rng.multivariate_normal(mean=(qx0, qy0, qz0, en0), cov=cov, size=sampling_pts)
        for q in pts:
            rez = resolution_fnc(q - np.array((qx0, qy0, qz0, en0)), rez_mat, r0)
            rez_conv += rez * signal(*q)
        rez_conv /= sampling_pts
        result.append(rez_conv)
    return np.squeeze(np.reshape(result, sz)) * N**2


# plot signal
N = 100
q1 = np.linspace(-1, 1, N)
en = np.linspace(-5, 5, N)


fig, axes = plt.subplots(ncols=2, nrows=2, sharex=True)
ax0 = axes[0, 0]
Q, E = np.meshgrid(q1, en)
img0 = ax0.pcolormesh(Q, E, signal(q1=q1, q2=0, q3=0, en=en).T, cmap="turbo", vmin=0, vmax=2000)
ax0.set_title("Bragg Peak")
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

ax2 = axes[0, 1]
ax2.plot(q1, signal(q1=q1, q2=0, q3=0, en=0), "-o")
# ax2.set_xlabel("Q")
ax2.set_ylabel("Intensity")
ax2.grid(alpha=0.6)
ax2.set_title("Bragg Peak E = 0")

# plot measurement
N = 40
q1 = np.linspace(-1, 1, N)
en = np.linspace(-5, 5, N)

ax1 = axes[1, 0]
Q, E = np.meshgrid(q1, en)
img1 = ax1.pcolormesh(Q, E, rez_conv_4d(q1=q1, q2=0, q3=0, en=en).T, cmap="turbo", vmin=0)
ax1.set_title("resolution convoluted")
ax1.grid(alpha=0.6)
ax1.set_xlabel("Q")
ax1.set_ylabel("E")
fig.colorbar(img1, ax=ax1)

ax3 = axes[1, 1]
ax3.plot(q1, rez_conv_4d(q1=q1, q2=0, q3=0, en=0), "-o")
ax3.set_xlabel("Q")
ax3.set_ylabel("Intensity")
ax3.grid(alpha=0.6)
ax3.set_title("resolution convoluted E = 0")

plt.tight_layout()
plt.show()
