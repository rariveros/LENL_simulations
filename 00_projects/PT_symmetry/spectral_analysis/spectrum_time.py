from back_process import *
from jacobians import *

def gaussian(x, A, x_0, sigma):
    return A * np.exp(- (x - x_0) ** 2 / (2 * sigma ** 2))

if __name__ == '__main__':
    eq = 'pdnlS'
    directory = 'C:/mnustes_science/simulation_data/FD/PT_extended/variando_distancia_rabi/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.0180/sigma=6.000/gamma=0.1850/dist=30.2930'

    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    x_grid = np.loadtxt(directory + '/X.txt', delimiter=',')
    t_grid = np.loadtxt(directory + '/T.txt', delimiter=',')
    params = np.loadtxt(directory + '/parameters.txt', delimiter=',')
    # [alpha, beta, gamma_0, dist,  mu, nu, sigma]
    dist = params[3]
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
    Nt = len(t_grid)
    each_T = 5
    ti = int(0.1 * Nt)

    Z_r = Z_r[ti::each_T, :]
    Z_i = Z_i[ti::each_T, :]
    t_grid = t_grid[ti::each_T]

    Nt = len(t_grid)
    colors = np.arange(0, 1, 1 / (Nt + 1))
    values_r = []
    values_i = []
    for i in range(Nt):
        print(str(i) + "/" + str(Nt))
        J = jacobians_FD(eq, [Z_r[i, :], Z_i[i, :]], [0], x_grid, [0], parameters, [DD])
        eigenvalues, eigenvectors = np.linalg.eig(J)
        #print(np.real(eigenvalues))
        arg_max = np.argmax(np.real(eigenvalues))
        x_array = np.real(eigenvalues) > -0.095
        values_r.append(np.real(eigenvalues))
        values_i.append(np.imag(eigenvalues))
        plt.scatter(np.real(eigenvalues[x_array]), np.imag(eigenvalues[x_array]), marker="D",c=[[colors[i], 0, 1 - colors[i]]], zorder=10)
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

    values_r = np.array(values_r)
    values_i = np.array(values_i)

    np.savetxt(directory + '/eigenvalues_r.txt', values_r, delimiter=',')
    np.savetxt(directory + '/eigenvalues_i.txt', values_i, delimiter=',')