from back_process import *
from jacobians import *

def gaussian(x, A, x_0, sigma):
    return A * np.exp(- (x - x_0) ** 2 / (2 * sigma ** 2))

cmap_reversed = matplotlib.cm.get_cmap('autumn_r')

if __name__ == '__main__':
    eq = 'pdnlS'
    dists = np.append(np.arange(13.91, 19.91, 0.25), np.array([18.410]))
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
        gamma_real = -(gaussian(x_grid, gamma_0, -dist / 2, sigma) - gaussian(x_grid, gamma_0, dist / 2, sigma))
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

    mode = np.exp(fav_val_01 * 0.01) * fav_vec_01 + np.exp(fav_val_02 * 0.01) * fav_vec_02
    steady_r = steady_i = 0 * Z_r[-1, :]

    #plt.plot(x_grid, np.imag(mode[:Nx]), color="purple")
    #plt.plot(x_grid, np.imag(mode[Nx:]), color="green")
    #plt.show()
    #plt.close()

    #plt.plot(x_grid, np.imag(fav_vec_02[:Nx]), color="b", linestyle="--")
    #plt.plot(x_grid, np.imag(fav_vec_02[Nx:]), color="r", linestyle="--")
    #plt.plot(x_grid, fav_vec_03[:Nx], color="b")
    #plt.plot(x_grid, fav_vec_03[Nx:], color="r")
    #plt.plot(x_grid, fav_vec_04[:Nx], color="b", linestyle="--")
    #plt.plot(x_grid, fav_vec_04[Nx:], color="r", linestyle="--")

    mode_evol = np.sqrt(steady_r ** 2 + steady_i ** 2 + 2 * mode[:Nx] * steady_r + 2 * steady_i * np.real(mode[Nx:]) + np.real(mode[:Nx]) ** 2 + np.real(mode[Nx:]) ** 2)
    mode_0 = np.sqrt(steady_r ** 2 + steady_i ** 2)
    eig_r = np.real(EIG)
    eig_i = np.imag(EIG)

    np.savetxt('spectral_data/mod_ev_zero.txt', np.abs(mode_evol), delimiter=',')
    np.savetxt('spectral_data/mod_0_zero.txt', np.abs(mode_0), delimiter=',')
    np.savetxt('spectral_data/eig_r_zero.txt', eig_r, delimiter=',')
    np.savetxt('spectral_data/eig_i_zero.txt', eig_i, delimiter=',')
    np.savetxt('spectral_data/DIST_zero.txt', DIST, delimiter=',')
    np.savetxt('spectral_data/x_grid_zero.txt', x_grid, delimiter=',')

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x_grid, np.sqrt(steady_r ** 2 + steady_i ** 2 + 2 * mode[:Nx] * steady_r + 2 * steady_i * np.real(mode[Nx:]) + np.real(mode[:Nx]) ** 2 + np.real(mode[Nx:]) ** 2), color="k", label="$|A_0 + e^{\lambda_c dt}A_c|$", linewidth=2)
    ax.plot(x_grid, np.sqrt(steady_r ** 2 + steady_i ** 2), color="k", label="$|A_0|$", linestyle=":", linewidth=2)
    #ax.set_aspect(0.6 * 300)
    ax.set_xlim(-75, 75)
    ax.grid(alpha=0.3)
    plt.xlabel("$x\ \\textrm{(mm)}$", fontsize=18)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)
    plt.savefig("plots/zero_mod_mode.png", dpi=200)
    plt.close()

    sc = plt.scatter(np.real(EIG) * 100, np.imag(EIG) * 100, edgecolors="0.2", c=DIST, zorder=10, cmap=np.flip(parula_map), s=60)
    plt.scatter(real_c * 100, imag_c * 100, edgecolors="0.2", c="r", zorder=20, s=60)
    cbar = plt.colorbar(sc)
    cbar.set_label('$d\ \\textrm{(mm)}$', rotation=0, size=22, labelpad=-27, y=1.1)
    cbar.ax.tick_params(labelsize=14)
    plt.grid(alpha=0.4, zorder=0)
    plt.vlines(0, -20, 20, colors="k", alpha=0.8)
    plt.hlines(0, -20, 20, colors="k", alpha=0.8)
    plt.xlim([-0.015 * 100, 0.003 * 100])
    plt.ylim([-0.1 * 100, 0.1 * 100])
    plt.xlabel("$\\textrm{Re}\{\lambda\}\ \\times 10^{-2}$", fontsize=18)
    plt.ylabel("$\\textrm{Im}\{\lambda\}\ \\times 10^{-2}$", fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    #plt.show()
    plt.savefig("plots/spectrum_zero.png", dpi=200)
    plt.close()