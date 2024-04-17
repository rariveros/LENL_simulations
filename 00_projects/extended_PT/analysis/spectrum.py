from back_process import *
from jacobians import *

if __name__ == '__main__':
    disco = 'D:/'
    eq = 'pdnlS'

    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FD'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    x_grid = np.loadtxt(directory + '/X.txt', delimiter=',')
    t_grid = np.loadtxt(directory + '/T.txt', delimiter=',')
    [alpha, beta, gamma_0, mu, nu, sigma] = np.loadtxt(directory + '/parameters.txt', delimiter=',')
    print("[alpha, beta, gamma_0, mu, nu, sigma]")
    print([alpha, beta, gamma_0, mu, nu, sigma])
    Nx = len(x_grid)
    dx = x_grid[1] - x_grid[0]
    phi = np.pi
    gamma_l = 1
    gamma_r = 1
    gamma_real = gamma_0
    gamma_img = 0
    gamma = [gamma_real, gamma_img]
    parameters = [alpha, beta, gamma, mu, nu]
    DD = sparse_DD_neumann(Nx, dx)
    for i in range(1, 2000, 100):
        J = jacobians_FD(eq, [Z_r[-i, :], Z_i[-i, :]], t_grid, x_grid, [0], parameters, [DD])
        plt.plot(Z_r[-i, :])
        plt.show()
        plt.close()
        eigenvalues, eigenvectors = np.linalg.eig(J)
        for j in range(len(eigenvectors)):
            #if np.real(eigenvalues[j]) > 0 and np.imag(eigenvalues[j]) != 0:
            print("Eigenvalue: " + str(eigenvalues[j]))
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