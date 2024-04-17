from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    x_grid = np.arange(-300, 300, 1)
    sigma = 45
    gaussian_L = np.exp(-(x_grid - 150) ** 2 / (2 * sigma ** 2))
    gaussian_R = np.exp(-(x_grid + 150) ** 2 / (2 * sigma ** 2))
    antigaussian_L = np.sign(x_grid - 150) * np.exp(-(x_grid - 150) ** 2 / (2 * sigma ** 2))
    antigaussian_R = np.sign(x_grid + 150) * np.exp(-(x_grid + 150) ** 2 / (2 * sigma ** 2))
    k_1 = 0.3
    k_2 = 0.3

    Z_r = antigaussian_L * np.sin(k_1 * x_grid) + gaussian_R * np.cos(k_1 * x_grid)
    Z_i = gaussian_L * np.cos(k_2 * x_grid) + antigaussian_R * np.sin(k_2 * x_grid)
    plt.plot(x_grid, Z_r, color="b")
    plt.plot(x_grid, Z_i, color="r")
    plt.plot(x_grid, np.sqrt(Z_r ** 2 + Z_i ** 2), color="k")
    plt.show()