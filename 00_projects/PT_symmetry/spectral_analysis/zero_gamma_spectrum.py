from back_process import *
from jacobians import *

def gaussian(x, A, x_0, sigma):
    return A * np.exp(- (x - x_0) ** 2 / (2 * sigma ** 2))

cmap_reversed = matplotlib.cm.get_cmap('autumn_r')

if __name__ == '__main__':
    eq = 'pdnlS'
    gammas = np.arange(0.1, 1, 0.01)
    real_maxs = []
    imag_maxs = []
    DIST = []
    EIG = []
    distss = []
    directory = "C:/mnustes_science/simulation_data/FD/PT_extended/variando_distancia/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.0180/sigma=6.000/gamma=0.1850/dist=30.6000"
    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    x_grid = np.loadtxt(directory + '/X.txt', delimiter=',')
    t_grid = np.loadtxt(directory + '/T.txt', delimiter=',')
    [alpha, beta, gamma_0_fake, dist, mu_0, nu, sigma] = np.loadtxt(directory + '/parameters.txt', delimiter=',')
    #nu = 0.030
    #alpha = 1
    #gamma_0 = 0.24
    #nu = 0.1
    #sigma = 16
    dist = 24
    for i in range(len(gammas)):
        print(" #################   dist=" + str(gammas[i]) + "   ###########")
        gamma_0 = gammas[i]
        Nx = len(x_grid)
        dx = x_grid[1] - x_grid[0]
        mu = mu_0 * np.ones(Nx)
        gamma_l = 1
        gamma_r = 1
        gamma_img = 0
        gamma_real = gamma_0 * np.sign(x_grid) #-(gaussian(x_grid, gamma_0, -dist / 2, sigma) - gaussian(x_grid, gamma_0, dist / 2, sigma))
        gamma = [gamma_real, gamma_img]
        parameters = [alpha, beta, gamma, mu, nu]
        DD = sparse_DD_neumann(Nx, dx)
        J = jacobians_FD(eq, [0 * Z_r[-1, :], 0 * Z_i[-1, :]], [0], x_grid, [0], parameters, [DD])
        eigenvalues, eigenvectors = np.linalg.eig(J)
        arg_max = np.argmax(np.real(eigenvalues))
        real_maxs.append(np.real(eigenvalues[arg_max]))
        imag_maxs.append(np.imag(eigenvalues[arg_max]))
        plt.scatter(gamma_0, np.imag(eigenvalues[arg_max]), color="r", edgecolors="k", zorder=10)
        plt.scatter(gamma_0, np.real(eigenvalues[arg_max]), color="b", edgecolors="k", zorder=10)
    plt.scatter(gamma_0, np.imag(eigenvalues[arg_max]), color="r", edgecolors="k", label="$\\textrm{Im}\{\lambda\}$")
    plt.scatter(gamma_0, np.real(eigenvalues[arg_max]), color="b", edgecolors="k", label="$\\textrm{Re}\{\lambda\}$")
    #plt.hlines(0, 150, 0, colors="k", alpha=0.8)
    #plt.vlines(30.294, -0.2, 0.2, colors="k", alpha=0.8, linestyles="--")
    #plt.xlim(0, 150)
    #plt.ylim(-0.15, 0.15)
    plt.legend(fontsize=15)
    plt.grid(alpha=0.4, zorder=0)
    plt.xlabel("$\gamma_0$", fontsize=20)
    plt.ylabel("$\lambda$", fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tight_layout()
    #plt.savefig("quiscient_drift.png", dpi=200)
    plt.show()
    plt.close()