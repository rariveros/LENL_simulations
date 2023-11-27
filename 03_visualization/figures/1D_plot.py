from back_process import *

if __name__ == '__main__':
    x_grid = np.arange(-50, 50, 0.25)
    Nx = len(x_grid)
    sigma = 15.3
    dist = 37
    gamma_0 = 1
    gamma_l = gamma_r = 1
    phi = np.pi
    y_grid = gamma_0 * (gamma_l * np.exp(- (x_grid - dist / 2) ** 2 / (2 * sigma ** 2)) + gamma_r * np.cos(phi) * np.exp(- (x_grid + dist / 2) ** 2 / (2 * sigma ** 2))) + 0.1 * (np.random.rand(Nx) - 0.5)

    matplotlib.rc('xtick', labelsize=15)
    matplotlib.rc('ytick', labelsize=15)

    fig, ax1 = plt.subplots()
    ax1.plot(x_grid, np.abs(y_grid), c="k", linewidth=2)
    ax1.set_xlabel("$x\ \\textrm{(mm)}$", fontsize=22)
    ax1.set_ylabel("$\langle \eta \\rangle\ \\textrm{(A.U.)}$", fontsize=22)
    ax1.set_xlim(x_grid[0], x_grid[-1])
    ax1.set_ylim(0, 1.2)#1.1 * np.amax(Z_modulo[I_R[0] + 10, :]))
    ax1.set_aspect(0.4 * (x_grid[-1] - x_grid[0])/1.1)
    ax1.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("response.png", dpi=200)
    plt.show()
    plt.close()