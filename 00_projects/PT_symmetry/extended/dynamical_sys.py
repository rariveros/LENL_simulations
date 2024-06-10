import matplotlib.pyplot as plt

from back_process import *
from jacobians import *

def gaussian(x, A, x_0, sigma):
    return A * np.exp(- (x - x_0) ** 2 / (2 * sigma ** 2))

if __name__ == '__main__':
    eq = 'pdnlS'
    disco = 'D:/'
    real_maxs = []
    imag_maxs = []
    DIST = []
    EIG = []
    parent_directory ="C:/mnustes_science/simulation_data/FD/PT_extended/variando_distancia_rabi/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.0180/sigma=6.000/gamma=0.1850"
    dist_dir = [name for name in os.listdir(parent_directory) if os.path.isdir(os.path.join(parent_directory, name))]
    j = 0
    colors = np.arange(0, 1, 1 / (len(dist_dir) + 1))
    for i in range(len(dist_dir)):
        print(" #################   " + dist_dir[i] + "    ###########")
        directory = parent_directory + "/" + dist_dir[i]
        [alpha, beta, gamma_0, dist, mu_0, nu, sigma] = np.loadtxt(directory + '/parameters.txt', delimiter=',')
        n = 1 - (4 / 30)
        #if dist % 2 == 0:
        beta = 0.004811649356064012
        if dist == 21:
            Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
            Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
            x_grid = np.loadtxt(directory + '/X.txt', delimiter=',')
            t_grid = np.loadtxt(directory + '/T.txt', delimiter=',')

            Z_r = filtro_superficie(Z_r, 5, "XY")
            Z_i = filtro_superficie(Z_i, 5, "XY")
            Nx = len(x_grid)
            Nt = len(t_grid)
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

            ZL_r = Z_r[:, :int(Nx / 2)]
            ZR_r = Z_r[:, int(Nx / 2) + 1:]
            ZL_i = Z_i[:, :int(Nx / 2)]
            ZR_i = Z_i[:, int(Nx / 2) + 1:]
            x_grid_R = x_grid[int(Nx / 2) + 1:]
            x_grid_L = x_grid[:int(Nx / 2)]
            CR_r = []
            CL_r = []
            CR_i = []
            CL_i = []
            for i in range(Nt):
                cR_r = integrate.simpson(ZR_r[i, :], x_grid_R)
                cL_r = integrate.simpson(ZL_r[i, :], x_grid_L)
                cR_i = integrate.simpson(ZR_i[i, :], x_grid_R)
                cL_i = integrate.simpson(ZL_i[i, :], x_grid_L)
                CR_r.append(cR_r)
                CL_r.append(cL_r)
                CR_i.append(cR_i)
                CL_i.append(cL_i)
            CR_r = np.array(CR_r)
            CL_r = np.array(CL_r)
            CR_i = np.array(CR_i)
            CL_i = np.array(CL_i)

            fig, ((ax1, ax2)) = plt.subplots(nrows=2, ncols=1)
            ax1.hlines(0, 0, 400, alpha=0.5, colors="k")
            ax1.plot(t_grid[int(n * Nt):] - t_grid[int(n * Nt)], CR_r[int(n * Nt):], color="purple")
            ax1.plot(t_grid[int(n * Nt):] - t_grid[int(n * Nt)], CR_i[int(n * Nt):], color="green")
            ax1.plot(t_grid[int(n * Nt):] - t_grid[int(n * Nt)], np.sqrt(CR_r[int(n * Nt):] ** 2 + CR_i[int(n * Nt):] ** 2), color="k")
            ax1.set_xlim(0, 400)
            ax1.grid(alpha=0.2)

            ax2.hlines(0, 0, 400, alpha=0.5, colors="k")
            ax2.plot(t_grid[int(n * Nt):] - t_grid[int(n * Nt)], CL_r[int(n * Nt):], color="purple")
            ax2.plot(t_grid[int(n * Nt):] - t_grid[int(n * Nt)], CL_i[int(n * Nt):], color="green")
            ax2.plot(t_grid[int(n * Nt):] - t_grid[int(n * Nt)], np.sqrt(CL_r[int(n * Nt):] ** 2 + CL_i[int(n * Nt):] ** 2), color="k")
            ax2.set_xlim(0, 400)
            ax2.grid(alpha=0.2)
            plt.show()
            plt.close()

            plt.plot(CR_r[int(n * Nt):], CL_i[int(n * Nt):], color=(colors[j], 0, 1 - colors[j]), label="$d = " + str(dist) + "$")
        j = j + 1
    plt.legend()
    plt.xlabel("$\psi_{RR}$", fontsize=20)
    plt.ylabel("$\psi_{LI}$", fontsize=20)
    plt.grid(alpha=0.2)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig("dynamical.png")
    plt.close()