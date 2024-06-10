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
    #parent_directory = "C:/mnustes_science/simulation_data/FD/PT_drift/windows/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.1000/sigma=16.000/gamma=0.2400"

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
        if i == 0:
            arg_max = np.argmax(np.real(eigenvalues))
            fav_val = eigenvalues[arg_max]
            fav_vec = eigenvectors[:, arg_max]
        elif i != 0:
            vars = []
            for j in range(len(eigenvalues)):
                var = np.abs(fav_val - eigenvalues[j])
                vars.append(var)
            arg_max = np.argmin(vars)
            fav_val = eigenvalues[arg_max]
        x_array = np.real(eigenvalues) > -0.095
        arg_otro = np.abs(np.imag(eigenvalues)) == np.abs(np.imag(fav_val))
        if dist == 30.295:
            fav_vec = eigenvectors[:, arg_otro]
            print(eigenvalues[arg_otro])
            fav_vec_01 = fav_vec[:, 0]
            #fav_vec_02 = fav_vec[:, 1]
            #fav_vec_03 = fav_vec[:, 2]
            fav_vec_02 = fav_vec[:, 3]
        real_maxs.append(np.real(eigenvalues[arg_max]))
        imag_maxs.append(np.imag(eigenvalues[arg_max]))
        EIG.append(eigenvalues)
        DIST.append([dist] * len(eigenvalues))
        #plt.scatter(np.real(eigenvalues[x_array]), np.imag(eigenvalues[x_array]), marker="D", c=[[colors[i], 0, 1 - colors[i]]], zorder=10, label="$" + dists[i] + "$")
        #plt.scatter(dist, np.imag(eigenvalues[x_array]), color="r", edgecolors="k", zorder=10)
        #plt.scatter(dist, np.real(eigenvalues[arg_max]), color="b", edgecolors="k", zorder=10)
        plt.scatter([dist] * len(eigenvalues[x_array]), np.imag(eigenvalues[x_array]), color="r", edgecolors="k", zorder=10)
        plt.scatter([dist] * len(eigenvalues[x_array]), np.real(eigenvalues[x_array]), color="b", edgecolors="k", zorder=10)
    plt.scatter(dist, np.imag(eigenvalues[arg_max]), color="r", edgecolors="k", label="$\\textrm{Im}\{\lambda\}$", zorder=10)
    plt.scatter(dist, np.real(eigenvalues[arg_max]), color="b", edgecolors="k", label="$\\textrm{Re}\{\lambda\}$", zorder=10)
    plt.vlines(30.294, -0.2, 0.2, colors="k", alpha=0.8, linestyles="--")
    plt.hlines(0, 80, 0, colors="k", alpha=0.8)
    plt.grid(alpha=0.4, zorder=0)
    plt.xlim(30, 150)
    plt.ylim(-0.1, 0.1)
    plt.legend(fontsize=15)
    plt.grid(alpha=0.4, zorder=0)
    plt.xlabel("$d\ \\textrm{(mm)}$", fontsize=20)
    plt.ylabel("$\lambda$", fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tight_layout()
    plt.savefig("stationary_draft.png", dpi=200)
    plt.close()

    plt.plot(x_grid, fav_vec_01[:Nx], color="b")
    plt.plot(x_grid, fav_vec_01[Nx:], color="r")
    plt.plot(x_grid, fav_vec_02[:Nx], color="b", linestyle="--")
    plt.plot(x_grid, fav_vec_02[Nx:], color="r", linestyle="--")
    plt.show()
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
    plt.savefig("spectrum.png", dpi=200)
    plt.close()

    X = np.array(DIST)[:, 0] < 40
    EIG = np.array(EIG)[X, :]
    sc = plt.scatter(np.real(EIG), np.imag(EIG), edgecolors="0.2", c=np.array(DIST)[X, :], zorder=10, cmap=parula_map)
    cbar = plt.colorbar(sc)
    cbar.set_label('$d\ \\textrm{(mm)}$', rotation=0, size=18, labelpad=-27, y=1.1)
    plt.grid(alpha=0.4, zorder=0)
    plt.hlines(-0.1, 0, 0, colors="k", alpha=0.8)
    plt.vlines(-0.2, 0.2, 0, colors="k", alpha=0.8)
    plt.xlim([-0.095, 0.01])
    plt.ylim([-0.15, 0.15])
    plt.xlabel("$\\textrm{Re}\{\lambda_i\}$", fontsize=18)
    plt.ylabel("$\\textrm{Im}\{\lambda_i\}$", fontsize=18)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tight_layout()
    plt.savefig("spectrum_zoom.png", dpi=200)
    plt.close()