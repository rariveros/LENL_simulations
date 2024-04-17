from back_process import *
from jacobians import *

def gaussian(x, A, x_0, sigma):
    return A * np.exp(- (x - x_0) ** 2 / (2 * sigma ** 2))

if __name__ == '__main__':
    eq = 'pdnlS'
    #dists = ["20.000", "22.000", "24.000", "26.000", "28.000"]
    dist = "30.2950"
    directory = "C:/mnustes_science/simulation_data/FD/PT_extended/variando_distancia/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.0180/sigma=6.000/gamma=0.1850/dist=" + dist
    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    x_grid = np.loadtxt(directory + '/X.txt', delimiter=',')
    t_grid = np.loadtxt(directory + '/T.txt', delimiter=',')
    [alpha, beta, gamma_0, dist, mu_0, nu, sigma] = np.loadtxt(directory + '/parameters.txt', delimiter=',')
    Nx = len(x_grid)
    dx = x_grid[1] - x_grid[0]
    Nx = len(x_grid)
    phi = np.pi
    mu = mu_0 * np.ones(Nx)
    gamma_l = 1
    gamma_r = 1
    gamma_img = 0
    gamma_real = -(gaussian(x_grid, gamma_0, -dist / 2, sigma) - gaussian(x_grid, gamma_0, dist / 2, sigma))
    gamma = [gamma_real, gamma_img]
    parameters = [alpha, beta, gamma, mu, nu]
    DD = sparse_DD_neumann(Nx, dx)
    t_grid_reduced = t_grid[-500:-1:10]
    Z_r = Z_r[-500:-1:10, :]
    Z_i = Z_i[-500:-1:10, :]
    Nt = len(t_grid_reduced)
    colors = np.arange(0, 1, 1 / (Nt + 1))
    for i in range(Nt):
        print(str(i) + "/" + str(Nt))
        J = jacobians_FD(eq, [Z_r[i, :], Z_i[i, :]], [0], x_grid, [0], parameters, [DD])
        eigenvalues, eigenvectors = np.linalg.eig(J)
        arg_max = np.argmax(np.real(eigenvalues))
        x_array = np.real(eigenvalues) > -0.05
        plt.scatter(np.real(eigenvalues[x_array]), np.imag(eigenvalues[x_array]), marker="D",c=(colors[i], 0, 1 - colors[i]), zorder=10)
    plt.grid(alpha=0.4, zorder=0)
    # plt.xlim([-0.1, 0.1])
    # plt.ylim([-1.2, 1.2])
    plt.xlabel("$\\textrm{Re}\{\lambda_i\}$", fontsize=18)
    plt.ylabel("$\\textrm{Im}\{\lambda_i\}$", fontsize=18)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.vlines(0, -1.2, 1.2, colors="k")
    plt.hlines(0, -0.1, 0.1, colors="k")
    plt.tight_layout()
    plt.show()
    plt.close()