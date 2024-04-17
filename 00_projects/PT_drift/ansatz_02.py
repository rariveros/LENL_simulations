from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    x_grid = np.arange(-300, 300, 1)
    sigma = 20
    d = 30
    gaussian_L = np.exp(-(x_grid - d / 2) ** 2 / (2 * sigma ** 2))
    gaussian_R = np.exp(-(x_grid + d / 2) ** 2 / (2 * sigma ** 2))
    k = 0.3

    f = gaussian_L + gaussian_R
    g = gaussian_R - gaussian_L
    Z_r = f * np.cos(k * x_grid)
    Z_i = g * np.sin(k * x_grid)
    plt.plot(x_grid, Z_r, color="b")
    plt.plot(x_grid, Z_i, color="r")
    plt.plot(x_grid, np.sqrt(Z_r ** 2 + Z_i ** 2), color="k")
    plt.show()