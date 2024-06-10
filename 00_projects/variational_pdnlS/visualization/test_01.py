from functions import *
from back_process import *
from time_integrators import *


if __name__ == '__main__':
    x_grid = np.arange(-50, 50, 0.1)
    dx = x_grid[1] - x_grid[0]
    alpha = 26.096
    beta = 1
    sigma_i = 15
    nu = -0.1
    mu = 0.075
    gamma_0 = 0.2
    X = [0]
    gamma = gamma_0 * (np.exp(-x_grid ** 2 / (2 * sigma_i ** 2)))
    gammaR = gamma_0 * (np.exp(-X[0] ** 2 / (2 * sigma_i ** 2)))
    phaseR = 0.5 * np.arccos(mu / gammaR)
    deltaR = (- nu + np.sqrt(gammaR ** 2 - mu ** 2))
    AR_i = np.sqrt(2 * deltaR) * np.cos(phaseR)
    Z_r = (AR_i / np.cosh(np.sqrt(deltaR / alpha) * (x_grid - (X[0]))))
    ddZ_r = np.append(0, np.diff(np.append(np.diff(Z_r) / dx, 0)) / dx)
    Z_i = (- (gamma + mu) - np.sqrt((gamma + mu) ** 2 - 4 * beta * Z_r * (nu * Z_r + beta * Z_r ** 3 + alpha * ddZ_r))) / (2 * Z_r * beta)
    plt.plot(x_grid, Z_r, color="b")
    plt.plot(x_grid, Z_i, color="r")
    plt.show()