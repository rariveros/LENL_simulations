import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    a_grid = np.arange(-3, 3, 0.01, dtype=np.complex_)
    mu_grid = np.arange(-3, 3, 0.01, dtype=np.complex_)
    A, MU = np.meshgrid(a_grid, mu_grid)
    eigen_plus = 0.5 * (-(MU * (A ** 2 - 1)) + np.sqrt(MU ** 2 * (A ** 2 - 1) ** 2 - 4))
    eigen_minus = 0.5 * (-(MU * (A ** 2 - 1)) - np.sqrt(MU ** 2 * (A ** 2 - 1) ** 2 - 4))
    fig, axs = plt.subplots(2, 2)
    ax1 = axs[0, 0]
    ax2 = axs[0, 1]
    ax3 = axs[1, 0]
    ax4 = axs[1, 1]
    cp_01 = ax1.pcolormesh(np.real(A), np.real(MU), np.real(eigen_plus), cmap=parula_map, shading='auto')
    cp_02 = ax2.pcolormesh(np.real(A), np.real(MU), np.imag(eigen_plus), cmap=parula_map, shading='auto')
    cp_03 = ax3.pcolormesh(np.real(A), np.real(MU), np.real(eigen_minus), cmap=parula_map, shading='auto')
    cp_04 = ax4.pcolormesh(np.real(A), np.real(MU), np.imag(eigen_minus), cmap=parula_map, shading='auto')

    fig.colorbar(cp_01)  # Add a colorbar to a plot
    ax1.set_title('$\\textrm{Re}\ \lambda_{-}$')
    ax1.set_ylabel('$\mu$', fontsize=12)

    fig.colorbar(cp_02)  # Add a colorbar to a plot
    ax2.set_title('$\\textrm{Im}\ \lambda_{-}$')
    ax2.set_ylabel('$\mu$', fontsize=12)

    fig.colorbar(cp_03)  # Add a colorbar to a plot
    ax3.set_title('$\\textrm{Re}\ \lambda_{-}$')
    ax3.set_xlabel('$a$', fontsize=12)
    ax3.set_ylabel('$\mu$', fontsize=12)

    fig.colorbar(cp_04)  # Add a colorbar to a plot
    ax4.set_title('$\\textrm{Im}\ \lambda_{-}$')
    ax4.set_xlabel('$a$', fontsize=12)
    ax4.set_ylabel('$\mu$', fontsize=12)

    plt.tight_layout()
    plt.savefig("eigen.png", dpi=300)
