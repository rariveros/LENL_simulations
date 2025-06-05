from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    x_grid = np.arange(-150, 150, 1)
    mu = 0.1
    gamma = 0.19
    nu = 0.1
    theta = 0.5 * np.arccos(mu / gamma)
    delta = np.sqrt(gamma ** 2 - mu ** 2) - nu
    R = np.cos(- delta * x_grid * np.sqrt(6.5))
    psi = R * np.exp(1j * theta)
    plt.plot(x_grid, np.real(psi), color="b", lw=2)
    plt.plot(x_grid, np.imag(psi), color="r", lw=2)
    plt.xlabel('$x\ \\textrm{(mm)}$', size='25')
    plt.xlim([x_grid[0], x_grid[-1]])
    plt.grid(linestyle='--', alpha=0.5)
    plt.show()