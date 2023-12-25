import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use("TkAgg")


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
    """
    Builds a design matrix using the given inputs.

    Args:
        x1 (numpy.ndarray): A 1D array of shape (n_samples,) containing the first set of samples.
        x2 (numpy.ndarray): A 1D array of shape (n_samples,) containing the second set of samples.
        bandwidth (float): The bandwidth parameter for the kernel.

    Returns:
        numpy.ndarray: A 2D array of shape (n_samples, n_samples) containing the design matrix.
    """
    return np.exp(
        -np.sum((x1[:, None] - x2[None]) ** 2, axis=-1) / (2 * bandwidth**2)
    )


def optimize_param(design_mat, y, regularizer, lr):  # lr: learning rate
    theta = np.zeros(len(y))
    # 単位ベクトル
    I = np.ones(len(y))
    # 　劣勾配法
    for _ in range(10000):
        diff = I + design_mat.dot(theta)
        grad = design_mat.T.dot(y)
        # 劣勾配
        subgrad = np.zeros(len(y))
        for i in range(len(y)):
            if diff[i] > 0:
                subgrad[i] = grad[i]
        theta -= lr * (subgrad + regularizer * design_mat.T.dot(theta))
    return theta


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
    plt.savefig("IS2-homework1.png")


x, y = generate_data(sample_size=200)
design_mat = build_design_mat(x, x, bandwidth=1.0)
theta = optimize_param(design_mat, y, regularizer=1, lr=0.0001)

visualize(theta, x, y)
