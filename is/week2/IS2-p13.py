import numpy as np
import matplotlib

matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

np.random.seed(1)


def generate_data(sample_size):
    a = np.linspace(0, 4 * np.pi, num=sample_size // 2)
    x = np.concatenate(
        [
            np.stack([a * np.cos(a), a * np.sin(a)], axis=1),
            np.stack([(a + np.pi) * np.cos(a), (a + np.pi) * np.sin(a)], axis=1),
        ]
    )
    x += np.random.random(size=x.shape)
    y = np.concatenate([np.ones(sample_size // 2), -np.ones(sample_size // 2)])
    return x, y


def build_design_mat(x1, x2, bandwidth):
    return np.exp(
        -np.sum((x1[:, None] - x2[None]) ** 2, axis=-1) / (2 * bandwidth**2)
    )


def optimize_param(design_mat, y, regularizer):
    return np.linalg.solve(
        design_mat.T.dot(design_mat) + regularizer * np.identity(len(y)),
        design_mat.T.dot(y),
    )


def visualize(theta, x, y, grid_size=100, x_min=-16, x_max=16):
    grid = np.linspace(x_min, x_max, grid_size)
    X, Y = np.meshgrid(grid, grid)
    mesh_grid = np.stack([np.ravel(X), np.ravel(Y)], axis=1)
    design_mat = build_design_mat(x, mesh_grid, bandwidth=1.0)
    plt.clf()
    plt.figure(figsize=(6, 6))
    plt.xlim(x_min, x_max)
    plt.ylim(x_min, x_max)
    plt.contourf(
        X,
        Y,
        np.reshape(np.sign(design_mat.T.dot(theta)), (grid_size, grid_size)),
        alpha=0.4,
        cmap=plt.cm.coolwarm,
    )
    plt.scatter(x[y == 1][:, 0], x[y == 1][:, 1], marker="$O$", c="blue")
    plt.scatter(x[y == -1][:, 0], x[y == -1][:, 1], marker="x", c="red")
    plt.savefig("l2-p14.png")


x, y = generate_data(sample_size=200)
design_mat = build_design_mat(x, x, bandwidth=1.0)
theta = optimize_param(design_mat, y, regularizer=0.01)

visualize(theta, x, y)
