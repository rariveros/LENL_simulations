from back_process import *
from jacobians import *

def gaussian(x, A, x_0, sigma):
    return A * np.exp(- (x - x_0) ** 2 / (2 * sigma ** 2))

def tanh_monticulo(x, A, x_0, sigma, emp):
    return A * (np.tanh(emp * (x - x_0 - sigma/2) + 1) - (np.tanh(emp * (x - x_0 + sigma/2) + 1)))

cmap_reversed = matplotlib.cm.get_cmap('autumn_r')

if __name__ == '__main__':
    eq = 'pdnlS'
    dists = np.arange(0.08, 0.25, 0.01) #np.append(np.arange(13.91, 19.91, 0.25), np.array([18.410]))
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
    [alpha, beta, gamma_0, dist, mu_0, nu, sigma] = np.loadtxt(directory + '/parameters.txt', delimiter=',')
    #nu = 0.030
    #alpha = 1
    #gamma_0 = 0.24
    #nu = 0.1
    #sigma = 16

    for i in range(len(dists)):
        print(" #################   dist=" + str(dists[i]) + "   ###########")
        dist = dists[i]
        Nx = len(x_grid)
        dx = x_grid[1] - x_grid[0]
        mu = mu_0 * np.ones(Nx)
        gamma_l = 1
        gamma_r = 1
        gamma_img = 0
        gamma_0 = dists[i]
        gamma_real = tanh_monticulo(x_grid, gamma_0, -10, sigma, 100) - tanh_monticulo(x_grid, gamma_0, 10, sigma, 100)
        #plt.plot(x_grid, gamma_real)
        #plt.show()
        gamma = [gamma_real, gamma_img]
        parameters = [alpha, beta, gamma, mu, nu]
        DD = sparse_DD_neumann(Nx, dx)
        J = jacobians_FD(eq, [0 * Z_r[-1, :], 0 * Z_i[-1, :]], [0], x_grid, [0], parameters, [DD])
        eigenvalues, eigenvectors = np.linalg.eig(J)
        arg_max = np.argmax(np.real(eigenvalues))
        x_array = np.real(eigenvalues) > -0.095
        PHI_R = []
        PHI_I = []
        EIG.append(eigenvalues)
        DIST.append([dist] * len(eigenvalues))
        distss.append(dist)
        #print(distss)
        real_maxs.append(np.real(eigenvalues[arg_max]))
        imag_maxs.append(np.imag(eigenvalues[arg_max]))
        plt.scatter(dist, np.imag(eigenvalues[arg_max]), color="r", edgecolors="k", zorder=10)
        plt.scatter(dist, np.real(eigenvalues[arg_max]), color="b", edgecolors="k", zorder=10)
        if dist == 18.410:
            real_c = np.real(eigenvalues)
            imag_c = np.imag(eigenvalues)
            fav_vec = eigenvectors[:, x_array]
            fav_val = eigenvalues[x_array]
            print(eigenvalues[x_array])
            fav_vec_01 = fav_vec[:, 0]
            fav_vec_02 = fav_vec[:, 1]
            fav_val_01 = fav_val[0]
            fav_val_02 = fav_val[1]
    #print(np.polyfit(distss, imag_maxs, 1))
    plt.scatter(dist, np.imag(eigenvalues[arg_max]), color="r", edgecolors="k", label="$\\textrm{Im}\{\lambda\}$")
    plt.scatter(dist, np.real(eigenvalues[arg_max]), color="b", edgecolors="k", label="$\\textrm{Re}\{\lambda\}$")
    plt.hlines(0, 150, 0, colors="k", alpha=0.8)
    #plt.vlines(30.294, -0.2, 0.2, colors="k", alpha=0.8, linestyles="--")
    plt.xlim(0, 150)
    #plt.ylim(-0.15, 0.15)
    plt.legend(fontsize=15)
    plt.grid(alpha=0.4, zorder=0)
    plt.xlabel("$d\ \\textrm{(mm)}$", fontsize=20)
    plt.ylabel("$\\textrm{ln} \lambda$", fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tight_layout()
    #plt.savefig("quiscient_drift.png", dpi=200)
    #plt.show()
    plt.close()

    sc = plt.scatter(np.real(EIG) * 100, np.imag(EIG) * 100, edgecolors="0.2", c=DIST, zorder=10, cmap=np.flip(parula_map), s=60)
    cbar = plt.colorbar(sc)
    cbar.set_label('$d\ \\textrm{(mm)}$', rotation=0, size=22, labelpad=-27, y=1.1)
    cbar.ax.tick_params(labelsize=14)
    plt.grid(alpha=0.4, zorder=0)
    plt.vlines(0, -20, 20, colors="k", alpha=0.8)
    plt.hlines(0, -20, 20, colors="k", alpha=0.8)
    #plt.xlim([-0.015 * 100, 0.003 * 100])
    #plt.ylim([-0.1 * 100, 0.1 * 100])
    plt.xlabel("$\\textrm{Re}\{\lambda\}\ \\times 10^{-2}$", fontsize=18)
    plt.ylabel("$\\textrm{Im}\{\lambda\}\ \\times 10^{-2}$", fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.show()
    #plt.savefig("plots/spectrum_zero.png", dpi=200)
    plt.close()