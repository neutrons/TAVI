import math

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse


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


def signal_gaussian_2d(x, y):
    x0, y0 = 0, 0
    sigma_x, sigma_y = 0.1, 0.1
    amp = 5
    return (
        amp
        * np.exp(-((x - x0) ** 2) / 2 / sigma_x**2)
        * np.exp(-((y - y0) ** 2) / 2 / sigma_y**2)
        / sigma_x
        / sigma_y
        / (2 * np.pi)
    )


def rez_conv_2d(x, y):

    n_sample = 20000
    pts_norm = np.random.default_rng().multivariate_normal(mean=(0, 0), cov=np.diag((1, 1)), size=n_sample)

    sz = x.shape
    result = []

    for x0, y0 in zip(x.flat, y.flat):
        rez_mat = resolution_matrix(x0, y0)
        cov = np.linalg.inv(rez_mat)

        eigenvalues, eigenvectors = np.linalg.eig(cov)
        pts = pts_norm @ eigenvectors @ np.diag(np.sqrt(eigenvalues)) @ eigenvectors.T
        pts += (x0, y0)

        rez_conv = np.mean(signal_gaussian_2d(pts[:, 0], pts[:, 1]))
        result.append(rez_conv)

    return np.reshape(result, sz)


if __name__ == "__main__":

    xmin, xmax = -20, 20
    ymin, ymax = -20, 20

    # -----------------plot signal-----------------------

    nx, ny = 1000, 2000

    x = np.linspace(xmin, xmax, nx)
    y = np.linspace(ymin, ymax, ny)
    vx, vy = np.meshgrid(x, y)

    fig, (ax0, ax1) = plt.subplots(ncols=2, sharex=True, sharey=True)
    signal = signal_gaussian_2d(vx, vy)
    integ_inten = np.sum(signal) * (xmax - xmin) / nx * (ymax - ymin) / ny
    print(f"integrated intensity of signal = {integ_inten}")
    img0 = ax0.pcolormesh(vx, vy, signal, cmap="turbo", vmin=0)
    ax0.set_title(f"sginal, integ. inten = {integ_inten:.3f}")
    ax0.grid(alpha=0.6)
    fig.colorbar(img0, ax=ax0)

    # ------------------ plot measurement ------------------------
    xstep, ystep = 2, 1
    nx, ny = math.ceil((xmax - xmin) / xstep), math.ceil((ymax - ymin) / ystep)

    x = np.linspace(xmin, xmax, nx)
    y = np.linspace(ymin, ymax, ny)
    vx, vy = np.meshgrid(x, y)

    measurement = rez_conv_2d(vx, vy)
    integ_inten = np.sum(measurement * xstep * ystep)
    print(f"integrated intensity of measurement = {integ_inten:.3f}")

    img1 = ax1.pcolormesh(vx, vy, measurement, cmap="turbo", vmin=0)
    fig.colorbar(img1, ax=ax1)
    ax1.set_title(f"resolution convoluted, integ. inten = {integ_inten:.3f}")
    ax1.grid(alpha=0.6)

    for i in range(3):
        ax1.add_artist(
            Ellipse(
                xy=(0, 0),
                width=5 * 2 * (i + 1),
                height=3 * 2 * (i + 1),
                angle=30,
                edgecolor=f"C{i+1}",
                facecolor="none",
                label=f"{i+1}-sigma",
                zorder=i + 2,
            )
        )

    ax1.legend()

    plt.tight_layout()
    plt.show()
