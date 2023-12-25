from __future__ import division, print_function

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from sklearn.model_selection import KFold

matplotlib.use("TkAgg")

np.random.seed(0)  # set the random seed for reproducibility


def generate_sample(xmin, xmax, sample_size):
    """
    Generate a sample dataset.

    Args:
        xmin (float): The minimum value of x.
        xmax (float): The maximum value of x.
        sample_size (int): The number of samples to generate.

    Returns:
        x (ndarray): The x values of the sample.
        target (ndarray): The target values of the sample.
    """
    x = np.linspace(start=xmin, stop=xmax, num=sample_size)
    pix = np.pi * x
    target = np.sin(pix) / pix + 0.1 * x
    noise = 0.05 * np.random.normal(loc=0.0, scale=1.0, size=sample_size)
    return x, target + noise


def calc_design_matrix(x, c, h):
    """
    Calculate the design matrix.

    Args:
        x (ndarray): The input values.
        c (ndarray): The centers of the basis functions.
        h (float): The bandwidth parameter.

    Returns:
        ndarray: The design matrix.
    """
    return np.exp(-((x[None] - c[:, None]) ** 2) / (2 * h**2))


# create sample
sample_size = 50
xmin, xmax = -3, 3
x, y = generate_sample(xmin=xmin, xmax=xmax, sample_size=sample_size)

# k-fold cross-validation
best_l, best_h = None, None


# implement here

# k-fold cross-validation
k = 5
l_list = [2**i for i in range(-30, 1)]
h_list = [2**i for i in range(-5, 6)]

best_l, best_h = 0, 0
best_score = float("inf")

# loop over all combinations of l and h
for l in l_list:
    for h in h_list:
        score = 0
        kf = KFold(n_splits=k, shuffle=True, random_state=0)
        for train_index, test_index in kf.split(x):
            x_train, x_test = x[train_index], x[test_index]
            y_train, y_test = y[train_index], y[test_index]
            k_train = calc_design_matrix(x_train, x_train, h)
            theta = np.linalg.solve(
                k_train.T.dot(k_train) + l * np.identity(len(k_train)),
                k_train.T.dot(y_train[:, None]),
            )
            k_test = calc_design_matrix(x_test, x_train, h)
            y_pred = k_test.T.dot(theta)
            score += np.mean((y_test - y_pred.flatten()) ** 2)
        score /= k
        if score < best_score:
            best_score = score
            best_l, best_h = l, h

print("Best (l, h): ({}, {})".format(best_l, best_h))

#
# visualization
#

# calculate design matrix
k = calc_design_matrix(x, x, best_h)

# solve the least square problem
theta = np.linalg.solve(k.T.dot(k) + best_l * np.identity(len(k)), k.T.dot(y[:, None]))

# create data to visualize the prediction
X = np.linspace(start=xmin, stop=xmax, num=5000)
K = calc_design_matrix(x, X, best_h)
prediction = K.dot(theta)

# plot the results
plt.clf()
plt.scatter(x, y, c="green", marker="o")
plt.plot(X, prediction)
plt.savefig("IS1-homework1.png")
