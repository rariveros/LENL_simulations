from back_process import *
from jacobians import *

if __name__ == '__main__':
    disco = 'aD:/'
    eq = 'pdnlS'

    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FD'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    x_grid = np.loadtxt(directory + '/X.txt', delimiter=',')
    t_grid = np.loadtxt(directory + '/T.txt', delimiter=',')
    [alpha, beta, gamma_0, dist, sigma,  mu, nu] = np.loadtxt(directory + '/parameters.txt', delimiter=',')
    print("[alpha, beta, gamma_0, mu, nu, sigma]")
    print([alpha, beta, gamma_0, mu, nu, sigma])
    Nx = len(x_grid)
    dx = x_grid[1] - x_grid[0]
    phi = np.pi
    gamma_complex = gamma_0 * (np.exp(- (x_grid - dist / 2) ** 2 / (2 * sigma ** 2)) + np.exp(1j * phi) * np.exp(- (x_grid + dist / 2) ** 2 / (2 * sigma ** 2)))
    gamma_real = np.real(gamma_complex)
    gamma_img = np.imag(gamma_complex)
    gamma = [gamma_real, gamma_img]
    parameters = [alpha, beta, gamma, mu, nu]
    DD = sparse_DD_neumann(Nx, dx)

    J = jacobians_FD(eq, [Z_r[-1, :], Z_i[-1, :]], t_grid, x_grid, [0], parameters, [DD])

    eigenvalues, eigenvectors = np.linalg.eig(J)
    sorted_indices = np.argsort(np.real(eigenvalues))[::-1]
    top_eigenvalues = eigenvalues[sorted_indices[:2]]
    top_eigenvectors = eigenvectors[:, sorted_indices[:2]]
    
    print(top_eigenvalues)

    plt.plot(x_grid, top_eigenvectors[:Nx, 0], color="b", lw=1)
    plt.plot(x_grid, top_eigenvectors[Nx:, 0], color="r", lw=1)
    plt.plot(x_grid, top_eigenvectors[:Nx, 1], color="b", lw=1)
    plt.plot(x_grid, top_eigenvectors[Nx:, 1], color="r", lw=1)
    #plt.plot(x_grid, Z_r[-1, :], color="b", lw=2)
    #plt.plot(x_grid, Z_i[-1, :], color="r", lw=2)
    plt.plot(x_grid, np.append(np.diff(Z_r[-1, :]), 0) / dx, color="b", lw=2)
    plt.plot(x_grid, np.append(np.diff(Z_i[-1, :]), 0) / dx, color="r", lw=2)
    plt.show()
    plt.close()

    plt.scatter(np.real(eigenvalues), np.imag(eigenvalues), c="b", zorder=10)
    plt.grid(alpha=0.4, zorder=0)
    #plt.xlim([-0.1, 0.1])
    #plt.ylim([-1.2, 1.2])
    plt.xlabel("$\\textrm{Re}\{\lambda_i\}$", fontsize=18)
    plt.ylabel("$\\textrm{Im}\{\lambda_i\}$", fontsize=18)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.vlines(0, -1.2, 1.2, colors="k")
    #plt.hlines(0, -0.1, 0.1, colors="k")
    plt.tight_layout()
    plt.show()
    plt.close()