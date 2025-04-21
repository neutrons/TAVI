import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse


def rotation_matrix_2d(theta_deg):
    theta = np.radians(theta_deg)
    c = np.cos(theta)
    s = np.sin(theta)
    return np.array([[c, -s], [s, c]])


mean = (1, 2)
sigma1, sigma2 = 5, 3
angle = -30
cov = rotation_matrix_2d(angle).T @ np.array([[sigma1**2, 0], [0, sigma2**2]]) @ rotation_matrix_2d(angle)
print(cov)

rng = np.random.default_rng()
x = rng.multivariate_normal(mean=mean, cov=cov, size=1000)
print(x.shape)  # (number of points, dimension)
print(x.mean(axis=0))  # close to mean=(1,2)
print(np.cov(x.T))


fig, ax = plt.subplots()
ax.plot(x[:, 0], x[:, 1], ".", alpha=0.5)
for i in range(3):
    ax.add_artist(
        Ellipse(
            xy=mean,
            width=sigma1 * 2 * (i + 1),
            height=sigma2 * 2 * (i + 1),
            angle=-angle,
            edgecolor=f"C{i+1}",
            facecolor="none",
            label=f"{i+1}-sigma",
        )
    )
ax.axis("equal")
ax.grid()
ax.legend()
plt.show()
