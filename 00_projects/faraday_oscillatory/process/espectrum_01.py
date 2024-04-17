from back_process import *
from jacobians import *
from scipy.sparse.linalg import eigs
if __name__ == '__main__':
    disco = 'D:/'
    eq = 'pdnlS'

    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FD/PDNLS_oscillatory/alpha=5.721/beta=1.000/mu=0.10'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    x_grid = np.loadtxt(directory + '/X.txt', delimiter=',')
    t_grid = np.loadtxt(directory + '/T.txt', delimiter=',')
    [alpha, beta, gamma_0, mu_0, nu, sigma] = np.loadtxt(directory + '/parameters.txt', delimiter=',')
    Nx = len(x_grid)
    dx = x_grid[1] - x_grid[0]
    print(sigma)
    gamma_1 = gamma_0 * np.exp(- x_grid ** 2 / (2 * sigma ** 2))
    gamma_2 = 0
    gamma = [gamma_1, gamma_2]
    mu = mu_0 * np.ones(Nx)
    mu[0:10] = 10
    mu[-10:-1] = 10
    parameters = [alpha, beta, gamma, mu, nu]
    DD = sparse_DD_neumann(Nx, dx)
    J = jacobians_FD(eq, [Z_r[-1, :], Z_i[-1, :]], t_grid, x_grid, [0], parameters, [DD])
    len_J = len(J)
    J = sparse.csr_matrix(J)
    eigenvalues, eigenvectors = eigs(J, k=len_J-2)
    print(eigenvalues)
    plt.scatter(np.real(eigenvalues), np.imag(eigenvalues), marker="D", c="r", zorder=10)
    plt.grid(alpha=0.4, zorder=0)
    #plt.xlim([-0.1, 0.1])
    #plt.ylim([-1.2, 1.2])
    plt.xlabel("$\\textrm{Re}\{\lambda_i\}$", fontsize=18)
    plt.ylabel("$\\textrm{Im}\{\lambda_i\}$", fontsize=18)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.vlines(0, -1.2, 1.2, colors="k")
    plt.hlines(0, -0.1, 0.1, colors="k")
    plt.tight_layout()
    plt.show()
    plt.close()