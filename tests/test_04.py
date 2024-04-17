from functions import *
from back_process import *
from time_integrators import *


if __name__ == '__main__':
    x_grid = np.arange(-40, 40, 0.1)
    dx = x_grid[1] - x_grid[0]
    alpha = 1
    beta = 1
    sigma_i = 3
    nu = -0.05
    mu = 0.1
    gamma = 0.105 #* np.exp(- x_grid ** 2 / (2 * sigma_i ** 2))

    A = np.arange(0.35, 0.45, 0.0025)
    b = np.arange(3, 4, 0.01)
    F = []
    for i in range(len(A)):
        Fi = []
        for j in range(len(b)):
            A_i = A[i]
            b_i = b[j]
            phi_01_i = (A_i / np.cosh(x_grid / b_i)) * np.sin(1.95 * np.pi)
            phi_02_i = (A_i / np.cosh(x_grid / b_i)) * np.cos(1.95 * np.pi)
            Dphi_01_i = np.append(np.diff(phi_01_i) / (dx), 0)
            Dphi_02_i = np.append(np.diff(phi_02_i) / (dx), 0)
            f_ij = (-(alpha / 2) * (Dphi_01_i ** 2 + Dphi_02_i ** 2) + (beta / 4) * (phi_01_i ** 2 + phi_02_i ** 2) ** 2 + (nu / 2) * (phi_01_i ** 2 + phi_02_i ** 2) + gamma * phi_01_i * phi_02_i) / (2 * mu)
            F_ij = integrate.simpson(f_ij, x_grid)
            Fi.append(F_ij)
        F.append(Fi)
    F = np.array(F)

    X, Y = np.meshgrid(b, A)
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    surf = ax.plot_surface(X, Y, F, cmap=parula_map, linewidth=0, antialiased=False)
    ax.set_xlabel("$b$")
    ax.set_ylabel("$A$")
    plt.show()