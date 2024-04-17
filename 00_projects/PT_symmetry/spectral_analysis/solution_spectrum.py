from back_process import *
from jacobians import *

def gaussian(x, A, x_0, sigma):
    return A * np.exp(- (x - x_0) ** 2 / (2 * sigma ** 2))

if __name__ == '__main__':
    eq = 'pdnlS'
    #dists = ["20.000", "22.000", "24.000", "26.000", "28.000"]
    dists = ["30.0000", "30.1000", "30.2000", "30.3000", "30.4000", "30.5000"]#, "30.6000",  "30.7000", "30.8000", "30.9000", "31.000"]
    colors = np.arange(0, 1, 1 / (len(dists) + 1))
    print(colors)
    for i in range(len(dists)):
        print(" #################   dist=" + dists[i] + "   ###########")
        directory = "C:/mnustes_science/simulation_data/FD/PT_extended/variando_distancia/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.0180/sigma=6.000/gamma=0.1850/dist=" + dists[i]
        Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
        Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
        x_grid = np.loadtxt(directory + '/X.txt', delimiter=',')
        t_grid = np.loadtxt(directory + '/T.txt', delimiter=',')
        [alpha, beta, gamma_0, dist, mu_0, nu, sigma] = np.loadtxt(directory + '/parameters.txt', delimiter=',')
        #print("[alpha, beta, gamma_0, dist,  mu, nu, sigma]")
        #print([alpha, beta, gamma_0, dist, mu, nu, sigma])
        #print(dist)
        Nx = len(x_grid)
        dx = x_grid[1] - x_grid[0]
        Nx = len(x_grid)
        phi = np.pi
        mu = mu_0 * np.ones(Nx)
        gamma_l = 1
        gamma_r = 1
        gamma_img = 0
        gamma_real = -(gaussian(x_grid, gamma_0, -dist / 2, sigma) - gaussian(x_grid, gamma_0, dist / 2, sigma))
        #plt.plot(x_grid, gamma_real)
        #plt.show()
        #plt.close()
        gamma = [gamma_real, gamma_img]
        parameters = [alpha, beta, gamma, mu, nu]
        DD = sparse_DD_neumann(Nx, dx)

        J = jacobians_FD(eq, [Z_r[-1, :], Z_i[-1, :]], [0], x_grid, [0], parameters, [DD])
        #plt.imshow(J)
        #plt.show()
        #plt.close()
        eigenvalues, eigenvectors = np.linalg.eig(J)
        #plt.scatter(np.real(eigenvalues), np.imag(eigenvalues), marker="D", zorder=10)
        arg_max = np.argmax(np.real(eigenvalues))
        x_array = np.real(eigenvalues) > -0.15
        #print((colors[i], 0, 1 - colors[i]))
        #print(eigenvalues[x_array])
        plt.scatter(np.real(eigenvalues[x_array]), np.imag(eigenvalues[x_array]), marker="D", c=(colors[i], 0, 1 - colors[i]), zorder=10, label="$" + dists[i] + "$")
        #plt.scatter(np.real(eigenvalues[arg_max]), np.imag(eigenvalues[arg_max]), marker="D", c=(colors[i], 0, 1 - colors[i]), zorder=10, label="$" + dists[i] + "$")
        plt.legend()
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