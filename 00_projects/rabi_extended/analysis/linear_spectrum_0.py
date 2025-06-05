from back_process import *
from jacobians import *

if __name__ == '__main__':
    disco = 'D:/'
    eq = 'pdnlS'

    alpha = 2 * 6.524
    beta = 1
    nu = 0.03
    mu = 0.10
    gamma_0 = 0.157
    sigma = 12.67

    dx = 1
    [xmin, xmax, dx] = [-100, 100, dx]
    x_grid = np.arange(xmin, xmax, dx)

    Nx = len(x_grid)
    dx = x_grid[1] - x_grid[0]
    phi = np.pi
    DD = sparse_DD_neumann(Nx, dx)
    dists = np.arange(20, 70)
    Z_r = np.zeros(Nx)
    Z_i = np.zeros(Nx)
    for i in range(len(dists)):
        dist = dists[i]
        gamma_complex = gamma_0 * (np.exp(- (x_grid - dist / 2) ** 2 / (2 * sigma ** 2)) + np.exp(1j * phi) * np.exp(- (x_grid + dist / 2) ** 2 / (2 * sigma ** 2)))
        gamma_real = np.real(gamma_complex)
        gamma_img = np.imag(gamma_complex)
        gamma = [gamma_real, gamma_img]
        parameters = [alpha, beta, gamma, mu, nu]

        J = jacobians_FD(eq, [Z_r, Z_i], [0], x_grid, [0], parameters, [DD])

        eigenvalues, eigenvectors = np.linalg.eig(J)
        sorted_indices = np.argsort(np.real(eigenvalues))[::-1]
        top_eigenvalues = eigenvalues[sorted_indices[:2]]
        plt.scatter(len(top_eigenvalues) * [dist], np.real(top_eigenvalues), c="b", zorder=10)
        plt.scatter(len(top_eigenvalues) * [dist], np.imag(top_eigenvalues), c="r", zorder=10)
    plt.grid(alpha=0.4, zorder=0)
    #plt.xlim([-0.1, 0.1])
    plt.ylim([-0.1, 0.1])
    plt.xlabel("$d$", fontsize=18)
    plt.ylabel("$\\lambda$", fontsize=18)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.vlines(0, -1.2, 1.2, colors="k")
    #plt.hlines(0, -0.1, 0.1, colors="k")
    plt.tight_layout()
    plt.show()
    plt.close()