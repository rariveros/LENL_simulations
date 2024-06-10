from functions import *
from back_process import *
from time_integrators import *


if __name__ == '__main__':
    x_grid = np.arange(-150, 150, 0.1)
    sigma_01 = 5
    sigma_02 = 15
    gauss_01 = np.exp(- x_grid ** 2 / (2 * sigma_01 ** 2))
    gauss_02 = np.exp(- (x_grid - 30) ** 2 / (2 * sigma_02 ** 2))
    plt.plot(x_grid, gauss_01 * gauss_02, color="k")
    plt.plot(x_grid, gauss_01, color="b")
    plt.plot(x_grid, gauss_02, color="r")

    plt.show()
