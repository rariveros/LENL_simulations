from back_process import *
from jacobians import *

def gaussian(x, A, x_0, sigma):
    return A * np.exp(- (x - x_0) ** 2 / (2 * sigma ** 2))

if __name__ == '__main__':
    eq = 'pdnlS'
    disco = 'D:/'
    #dists = ["20.000", "22.000", "24.000", "26.000", "28.000"]
    real_maxs = []
    imag_maxs = []
    DIST = []
    EIG = []
    parent_directory ="C:/mnustes_science/simulation_data/FD/PT_extended/variando_distancia/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.0180/sigma=6.000/gamma=0.1850"
    dist_dir = [name for name in os.listdir(parent_directory) if os.path.isdir(os.path.join(parent_directory, name))]
    colors = np.arange(0, 1, 1 / (len(dist_dir) + 1))
    for i in range(len(dist_dir)):
        print(" #################   " + dist_dir[i] + "    ###########")
        directory = parent_directory + "/" + dist_dir[i]
        Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
        Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
        x_grid = np.loadtxt(directory + '/X.txt', delimiter=',')
        t_grid = np.loadtxt(directory + '/T.txt', delimiter=',')
        [alpha, beta, gamma_0, dist, mu_0, nu, sigma] = np.loadtxt(directory + '/parameters.txt', delimiter=',')

        #plt.plot(x_grid, Z_r[-1, :], label="R")
        #plt.plot(x_grid, Z_i[-1, :], label="I")
        #plt.plot(x_grid, np.sqrt(Z_i[-1, :] ** 2 + Z_i[-1, :] ** 2), color="k")
        #plt.legend()

        #plt.show()
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

        J = jacobians_FD(eq, [Z_r[-1, :], Z_i[-1, :]], [0], x_grid, [0], parameters, [DD])
        eigenvalues, eigenvectors = np.linalg.eig(J)
        arg_max = np.argmax(np.real(eigenvalues))
        x_array = np.real(eigenvalues) > -0.095
        real_maxs.append(np.real(eigenvalues[arg_max]))
        imag_maxs.append(np.imag(eigenvalues[arg_max]))
        EIG.append(eigenvalues)
        DIST.append([dist] * len(eigenvalues))
        if dist == 30.295:
            fav_vec = eigenvectors[:, arg_max]
            fav_val = eigenvalues[arg_max]
            steady_r = Z_r[-1, :]
            steady_i = Z_i[-1, :]
        #plt.scatter(np.real(eigenvalues[x_array]), np.imag(eigenvalues[x_array]), marker="D", c=[[colors[i], 0, 1 - colors[i]]], zorder=10, label="$" + dists[i] + "$")
        plt.scatter([dist] * len(eigenvalues[x_array]), np.imag(eigenvalues[x_array]), color="r", edgecolors="k", zorder=10)
        plt.scatter([dist] * len(eigenvalues[x_array]), np.real(eigenvalues[x_array]), color="b", edgecolors="k", zorder=10)
    EIG = np.array(EIG)
    plt.scatter([dist] * len(eigenvalues[x_array]), np.imag(eigenvalues[x_array]), color="r", edgecolors="k", label="$\\textrm{Im}\{\lambda\}$", zorder=10)
    plt.scatter([dist] * len(eigenvalues[x_array]), np.real(eigenvalues[x_array]), color="b", edgecolors="k", label="$\\textrm{Re}\{\lambda\}$", zorder=10)
    plt.vlines(30.294, -0.2, 0.2, colors="k", alpha=0.8, linestyles="--")
    plt.hlines(0, 80, 0, colors="k", alpha=0.8)
    plt.grid(alpha=0.4, zorder=0)
    plt.xlim(30, 31)
    plt.ylim(-0.15, 0.15)
    plt.legend(fontsize=15)
    plt.grid(alpha=0.4, zorder=0)
    plt.xlabel("$d\ \\textrm{(mm)}$", fontsize=20)
    plt.ylabel("$\lambda$", fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tight_layout()
    plt.savefig("plots/stationary_draft.png", dpi=200)
    plt.close()

    delta = np.exp(fav_val * 0.01)
    mode_evol = np.sqrt(steady_r ** 2 + steady_i ** 2 + 2 * delta * steady_r * np.real(fav_vec[:Nx]) + 2 * delta * steady_i * np.real(fav_vec[Nx:]) + delta ** 2 * np.real(fav_vec[:Nx]) ** 2 + delta ** 2 * np.real(fav_vec[Nx:]) ** 2)
    mode_0 = np.sqrt(steady_r ** 2 + steady_i ** 2)
    eig_r = np.real(EIG)
    eig_i = np.imag(EIG)

    np.savetxt('spectral_data/mod_ev_sol.txt', np.abs(mode_evol), delimiter=',')
    np.savetxt('spectral_data/mod_0_sol.txt', np.abs(mode_0), delimiter=',')
    np.savetxt('spectral_data/eig_r_sol.txt', eig_r, delimiter=',')
    np.savetxt('spectral_data/eig_i_sol.txt', eig_i, delimiter=',')
    np.savetxt('spectral_data/DIST_sol.txt', DIST, delimiter=',')
    np.savetxt('spectral_data/x_grid_sol.txt', x_grid, delimiter=',')

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x_grid, steady_r + delta * np.real(fav_vec[:Nx]), color="r", alpha=0.75,  label="$\\textrm{Re}(A_{s} + \delta A)$")
    ax.plot(x_grid, steady_i + delta * np.real(fav_vec[Nx:]), color="b", alpha=0.75, label="$\\textrm{Im}(A_{s} + \delta A)$")
    ax.plot(x_grid, steady_r, color="r")
    ax.plot(x_grid, steady_i, color="b")
    ax.set_aspect(0.6 * 300)
    ax.set_xlim(-75, 75)
    ax.grid(alpha=0.3)
    plt.xlabel("$x\ \\textrm{(mm)}$", fontsize=18)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig("plots/unstable_mode.png", dpi=300)
    plt.close()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x_grid, np.sqrt(steady_r ** 2 + steady_i ** 2 + 2 * delta * steady_r * np.real(fav_vec[:Nx]) + 2 * delta * steady_i * np.real(fav_vec[Nx:]) + delta ** 2 * np.real(fav_vec[:Nx]) ** 2 + delta ** 2 * np.real(fav_vec[Nx:]) ** 2), color="k", label="$|A_s + e^{\lambda_c dt}A_c|$", linewidth=2)
    ax.plot(x_grid, np.sqrt(steady_r ** 2 + steady_i ** 2), color="k", label="$|A_s|$", linestyle=":", linewidth=2)
    #ax.set_aspect(0.6 * 300)
    ax.set_xlim(-75, 75)
    ax.grid(alpha=0.3)
    plt.xlabel("$x\ \\textrm{(mm)}$", fontsize=18)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)
    plt.savefig("plots/unstable_mod_mode.png", dpi=200)
    plt.close()

    sc = plt.scatter(np.real(EIG), np.imag(EIG), edgecolors="0.2", c=DIST, zorder=10, cmap=parula_map)
    cbar = plt.colorbar(sc)
    cbar.set_label('$d\ \\textrm{(mm)}$', rotation=0, size=18, labelpad=-27, y=1.1)
    plt.grid(alpha=0.4, zorder=0)
    ## plt.xlim([-0.1, 0.1])
    ## plt.ylim([-1.2, 1.2])
    plt.xlabel("$\\textrm{Re}\{\lambda_i\}$", fontsize=18)
    plt.ylabel("$\\textrm{Im}\{\lambda_i\}$", fontsize=18)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tight_layout()
    plt.savefig("plots/spectrum_steady.png", dpi=200)
    plt.close()

    sc = plt.scatter(np.real(EIG) * 100, np.imag(EIG) * 100, edgecolors="0.2", c=DIST, zorder=10, cmap=parula_map, s=60)
    cbar = plt.colorbar(sc)
    cbar.set_label('$d\ \\textrm{(mm)}$', rotation=0, size=22, labelpad=-27, y=1.1)
    cbar.ax.tick_params(labelsize=14)
    plt.grid(alpha=0.4, zorder=0)
    plt.vlines(0, -20, 20, colors="k", alpha=0.8)
    plt.hlines(0, -20, 20, colors="k", alpha=0.8)
    plt.xlim([-0.08 * 100, 0.01 * 100])
    plt.ylim([-0.05 * 100, 0.05 * 100])
    plt.xlabel("$\\textrm{Re}\{\lambda\}\ \\times 10^{-2}$", fontsize=18)
    plt.ylabel("$\\textrm{Im}\{\lambda\}\ \\times 10^{-2}$", fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.savefig("plots/spectrum_steady_zoom.png", dpi=200)
    plt.close()