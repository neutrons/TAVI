import math
from time import time

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse


def timeit(func):
    def wrap_func(*args, **kwargs):
        t1 = time()
        result = func(*args, **kwargs)
        t2 = time()
        print(f"Function {func.__name__!r} executed in {(t2 - t1):.4f} s")
        return result

    return wrap_func


def gaussian_2d(cen, sigma, angle):
    "return center and covariance matrix of a 2D Gaussian distribution"

    def rotation_matrix_2d(theta_deg):
        theta = np.radians(theta_deg)
        c = np.cos(theta)
        s = np.sin(theta)
        return np.array([[c, -s], [s, c]])

    cov = rotation_matrix_2d(angle).T @ np.diag([sigma[0] ** 2, sigma[1] ** 2]) @ rotation_matrix_2d(angle)
    return (cen, cov)


@timeit
def gaussian_sampling(cen, cov, n_sample=10000):
    "directly generation using np.random"
    return np.random.default_rng().multivariate_normal(mean=cen, cov=cov, size=n_sample)


@timeit
def similarity_transformation(cen, cov, n_sample=10000):
    "using similarity transformation"
    x = np.random.default_rng().multivariate_normal(mean=(0, 0), cov=np.diag((1, 1)), size=n_sample)
    eigenvalues, eigenvectors = np.linalg.eig(cov)
    y = x @ eigenvectors @ np.diag(np.sqrt(eigenvalues)) @ eigenvectors.T
    y += cen
    return y


if __name__ == "__main__":
    cen = np.array([1, 2])
    sigma = np.array([5, 3])
    angle = -30

    cen, cov = gaussian_2d(cen, sigma, angle)
    print(f"2D Gaussian, center= {cen}, covariant matrix = ", *cov)

    n_sample = 10000
    x = gaussian_sampling(cen, cov, n_sample)
    # print(f"shape of x = {x.shape}")  # (number of points, dimension)
    print(f"mean of x = {x.mean(axis=0)}, cov of x = ", *np.cov(x.T))

    y = similarity_transformation(cen, cov, n_sample)
    print(f"mean of y = {y.mean(axis=0)}, cov of y = ", *np.cov(y.T))

    fig, axes = plt.subplots(ncols=2, nrows=2, sharex=True, sharey=True)

    # ------------------- fig 1------------------------
    ax0 = axes[0, 0]
    ax0.plot(x[:, 0], x[:, 1], ".", alpha=0.2)
    for i in range(3):
        ax0.add_artist(
            Ellipse(
                xy=cen,
                width=sigma[0] * 2 * (i + 1),
                height=sigma[1] * 2 * (i + 1),
                angle=-angle,
                edgecolor=f"C{i + 1}",
                facecolor="none",
                label=f"{i + 1}-sigma",
                zorder=i + 2,
            )
        )
    # ax0.axis("equal")
    ax0.set_xlabel("X")
    ax0.set_ylabel("Y")
    ax0.grid(alpha=0.6)
    ax0.set_title("Gaussian sampling")
    ax0.legend()

    # ------------------- fig 2 histo -------------------------
    xmin, xmax, xstep = -20, 22, 1
    ymin, ymax, ystep = -15, 17, 0.2
    xbins, ybins = math.ceil((xmax - xmin) / xstep), math.ceil((ymax - ymin) / ystep)
    xhist, (xedges, yedges) = np.histogramdd(x, bins=(xbins, ybins), range=((xmin, xmax), (ymin, ymax)), density=True)

    ax1 = axes[0, 1]
    im = ax1.pcolormesh(xedges, yedges, xhist.T, vmin=0, cmap="turbo", shading="flat")
    ax1.set_xlabel("X")
    ax1.set_ylabel("Y")
    ax1.grid(alpha=0.6)
    ax1.set_title("Gaussian sampling")
    fig.colorbar(im)

    # -------------------- fig 3 ------------------------
    ax2 = axes[1, 0]
    ax2.plot(y[:, 0], y[:, 1], ".", alpha=0.2)
    for i in range(3):
        ax2.add_artist(
            Ellipse(
                xy=cen,
                width=sigma[0] * 2 * (i + 1),
                height=sigma[1] * 2 * (i + 1),
                angle=-angle,
                edgecolor=f"C{i + 1}",
                facecolor="none",
                label=f"{i + 1}-sigma",
                zorder=i + 2,
            )
        )
    # ax0.axis("equal")
    ax2.set_xlabel("X")
    ax2.set_ylabel("Y")
    ax2.grid(alpha=0.6)
    ax2.set_title("Similarity transformation")
    ax2.legend()

    # ----------------------fig 4 histo----------------------
    yhist, (xedges, yedges) = np.histogramdd(y, bins=(xbins, ybins), range=((xmin, xmax), (ymin, ymax)), density=True)

    ax3 = axes[1, 1]
    im = ax3.pcolormesh(xedges, yedges, yhist.T, vmin=0, cmap="turbo")
    ax3.set_xlabel("X")
    ax3.set_ylabel("Y")
    ax3.grid(alpha=0.6)
    fig.colorbar(im)
    ax3.set_title("Similarity transformation")
    plt.tight_layout()

    # ---------------------- plot cuts ----------------
    fig, axes = plt.subplots(ncols=2, nrows=2, sharey="row", sharex="col")
    ax0 = axes[0, 0]
    ax0.plot(xedges[:-1], xhist[:, np.where(np.isclose(yedges, cen[1]))[0]], "o")

    sigma = np.sqrt(cov[0, 0])
    ax0.plot(
        xedges,
        1
        / (sigma * np.sqrt(2 * np.pi))
        * np.exp(-((xedges - cen[0]) ** 2) / (2 * sigma**2))
        / (np.sqrt(cov[1, 1]) * np.sqrt(2 * np.pi)),
        "-",
    )
    ax0.grid(alpha=0.6)
    ax0.set_xlabel("X")
    ax0.set_ylabel("counts")
    ax0.set_title("Gaussian sampling")

    # ------------------------------------
    ax1 = axes[0, 1]
    ax1.plot(yedges[:-1], xhist[np.where(np.isclose(xedges, cen[0]))[0], :].T, "o")
    sigma = np.sqrt(cov[1, 1])
    ax1.plot(
        yedges,
        1
        / (sigma * np.sqrt(2 * np.pi))
        * np.exp(-((yedges - cen[1]) ** 2) / (2 * sigma**2))
        / (np.sqrt(cov[0, 0]) * np.sqrt(2 * np.pi)),
        "-",
    )
    ax1.grid(alpha=0.6)
    ax1.set_xlabel("Y")
    ax1.set_ylabel("counts")
    ax1.set_title("Gaussian sampling")
    # ------------------------------------
    ax2 = axes[1, 0]
    ax2.plot(xedges[:-1], yhist[:, np.where(np.isclose(yedges, cen[1]))[0]], "o")

    sigma = np.sqrt(cov[0, 0])
    ax2.plot(
        xedges,
        1
        / (sigma * np.sqrt(2 * np.pi))
        * np.exp(-((xedges - cen[0]) ** 2) / (2 * sigma**2))
        / (np.sqrt(cov[1, 1]) * np.sqrt(2 * np.pi)),
        "-",
    )
    ax2.grid(alpha=0.6)
    ax2.set_xlabel("X")
    ax2.set_ylabel("counts")
    ax2.set_title("Similarity transformation")

    # ------------------------------------
    ax3 = axes[1, 1]
    ax3.plot(yedges[:-1], yhist[np.where(np.isclose(xedges, cen[0]))[0], :].T, "o")

    sigma = np.sqrt(cov[1, 1])
    ax3.plot(
        yedges,
        1
        / (sigma * np.sqrt(2 * np.pi))
        * np.exp(-((yedges - cen[1]) ** 2) / (2 * sigma**2))
        / (np.sqrt(cov[0, 0]) * np.sqrt(2 * np.pi)),
        "-",
    )
    ax3.grid(alpha=0.6)
    ax3.set_xlabel("Y")
    ax3.set_ylabel("counts")
    ax3.set_title("Similarity transformation")

    plt.tight_layout()
    plt.show()
