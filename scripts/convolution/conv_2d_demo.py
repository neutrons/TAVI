import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse


def signal_gaussian_2d(x, y):
    x0, y0 = 0, 0
    sigma_x, sigma_y = 0.1, 0.1
    return (
        np.exp(-((x - x0) ** 2) / 2 / sigma_x**2)
        * np.exp(-((y - y0) ** 2) / 2 / sigma_y**2)
        / sigma_x
        / sigma_y
        / (2 * np.pi)
    )


def resolution_matrix(x, y):

    def rotation_matrix_2d(theta_deg):
        theta = np.radians(theta_deg)
        c = np.cos(theta)
        s = np.sin(theta)
        return np.array([[c, -s], [s, c]])

    sigma1, sigma2 = 5, 3
    angle = -30
    rez_mat = (
        rotation_matrix_2d(angle).T @ np.array([[1 / sigma1**2, 0], [0, 1 / sigma2**2]]) @ rotation_matrix_2d(angle)
    )
    return rez_mat


def resolution_fnc(q, m):
    return np.sqrt(np.linalg.det(m)) * np.exp(-q.T @ m @ q / 2) / (2 * np.pi)


def rez_conv_2d(x, y):
    sz = x.shape
    result = []
    for x0, y0 in zip(x.flat, y.flat):
        rez_mat = resolution_matrix(x0, y0)
        cov = np.linalg.inv(rez_mat)

        sampling_pts = 2000
        rez_conv = 0
        rng = np.random.default_rng()
        pts = rng.multivariate_normal(mean=(x0, y0), cov=cov, size=sampling_pts)
        for q in pts:
            rez = resolution_fnc(q - np.array((x0, y0)), rez_mat)
            rez_conv += rez * signal_gaussian_2d(*q)
        rez_conv /= sampling_pts
        result.append(rez_conv)
    return np.reshape(result, sz) * Nx * Ny


# plot signal
Nx, Ny = 100, 200
X = np.linspace(-10, 10, Nx)
Y = np.linspace(-10, 10, Ny)
X, Y = np.meshgrid(X, Y)

fig, (ax0, ax1) = plt.subplots(ncols=2, sharex=True, sharey=True)
img0 = ax0.pcolormesh(X, Y, signal_gaussian_2d(X, Y), cmap="turbo", vmin=0, vmax=8)
ax0.set_title("sginal")
ax0.grid(alpha=0.6)
fig.colorbar(img0, ax=ax0)


# plot measurement
Nx, Ny = 10, 20
X = np.linspace(-10, 10, Nx)
Y = np.linspace(-10, 10, Ny)
X, Y = np.meshgrid(X, Y)
img1 = ax1.pcolormesh(
    X,
    Y,
    rez_conv_2d(X, Y),
    cmap="turbo",
    vmin=0,
)
ax1.set_title("resolution convoluted")
ax1.grid(alpha=0.6)
fig.colorbar(img1, ax=ax1)
ax1.add_artist(
    Ellipse(
        xy=(0, 0),
        width=5 * 2,
        height=3 * 2,
        angle=30,
        edgecolor="w",
        facecolor="none",
        label="1-sigma ellipse",
    )
)
ax1.legend()

plt.show()
