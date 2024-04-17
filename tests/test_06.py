from functions import *
from back_process import *
from time_integrators import *


if __name__ == '__main__':

    x_grid = np.arange(-25, 25, 0.025)
    dx = x_grid[1] - x_grid[0]
    a = np.arange(0.1, 1, 0.0025)
    alpha = 1
    V = +x_grid ** 2 / 2
    F = []
    for i in range(len(a)):
        phi_01_i = np.sqrt((16 * np.sqrt(a[i]) / (5 * np.pi))) / (1 + a[i] * x_grid ** 2) ** 2
        phi_02_i = np.sqrt((16 * np.sqrt(a[i]) / (5 * np.pi))) / (1 + a[i] * x_grid ** 2) ** 2

        Dphi_01_i = np.append(np.diff(phi_01_i) / (dx), 0) #np.append(np.diff(np.append(np.diff(phi_01_i) / (dx), 0)) / dx, 0)
        Dphi_02_i = np.append(np.diff(phi_02_i) / (dx), 0) #np.append(np.diff(np.append(np.diff(phi_01_i) / (dx), 0)) / dx, 0)

        f_ij = (alpha / 2) * (Dphi_01_i ** 2 + Dphi_02_i ** 2) + V * (phi_01_i ** 2 + phi_02_i ** 2)
        if a[i] == 1:
            plt.plot(x_grid, f_ij)
            plt.show()
            plt.close()
        F_i = integrate.simpson(f_ij, x_grid)
        F.append(F_i)
    F = np.array(F)

    plt.plot(a, F)
    plt.show()