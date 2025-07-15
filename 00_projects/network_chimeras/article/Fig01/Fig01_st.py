from back_process import *

if __name__ == '__main__':

    dir_01 = "D:/mnustes_science/simulation_data/FD/network_chimeras/FIG01/erdos_renyi/k=0.0150/mean_degree=18.00"
    U_light_01 = np.loadtxt(dir_01 + '/U_lyap.txt', delimiter=',')
    V_light_01 = np.loadtxt(dir_01 + '/V_lyap.txt', delimiter=',')

    t_light = np.loadtxt(dir_01 + '/T_lyap.txt', delimiter=',')
    x_grid = np.loadtxt(dir_01 + '/X.txt', delimiter=',')
    parameters_np = np.loadtxt(dir_01 + '/params.txt', delimiter=',')

    phase = np.arctan2(V_light_01, U_light_01)
    phase_light_wraped = phase
    phase_light = np.unwrap(phase_light_wraped, axis=0)
    module = np.sqrt(U_light_01 ** 2 + V_light_01 ** 2)

    dt = t_light[1] - t_light[0]
    Nt = len(t_light)
    Nx = len(x_grid)

    t_init = 4960
    t_final = 5000
    i_0 = np.argmin(np.abs(t_light - t_init))
    i_f = np.argmin(np.abs(t_light - t_final))
    power_threshold = 1.5

    file = dir

    arg_lyap_pos = np.loadtxt(file + '/arg_lyap.txt', delimiter=',')
    x = np.loadtxt(file + '/x_lyap.txt', delimiter=',')
    lyap_mean = np.loadtxt(file + '/lyap_spectrum.txt', delimiter=',')
    Q_pos = np.loadtxt(file + '/Q_pos.txt', delimiter=',')

    arg_lyap_pos = lyap_mean >= 0.0
    args_lyap_modes = np.argsort(Q_pos)

    #### SPATIOTEMPORAL DIAGRAMS ####
    N_interval = 2
    fig, ((ax01, ax02), (ax03, ax04)) = plt.subplots(2, 2, figsize=(10, 8))
    cax_01 = ax01.pcolormesh(x_grid, t_light[i_0:i_f:N_interval], U_light_01[i_0:i_f:N_interval , :], cmap="turbo", shading='auto')
    cbar = fig.colorbar(cax_01)
    cbar.ax.tick_params(labelsize=13)
    cbar.set_label('$u(t)$', rotation=0, size=20, labelpad=-50, y=1.1)
    ax01.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
    ax01.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=False, labeltop=False, labelbottom=False)
    ax01.set_ylabel("$t$", fontsize=20)

    cax_02 = ax02.pcolormesh(x_grid, t_light[i_0:i_f:N_interval ], module[i_0:i_f:N_interval , :], cmap="turbo", shading='auto')
    cbar = fig.colorbar(cax_02)
    cbar.ax.tick_params(labelsize=13)
    cbar.set_label('$r(t)$', rotation=0, size=20, labelpad=-50, y=1.1)
    ax02.tick_params(axis="y", direction="in", labelsize=15, left=False, right=True, labelleft=False, labelright=False)
    ax02.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=False, labeltop=False, labelbottom=False)

    cax_03 = ax03.pcolormesh(x_grid, t_light[i_0:i_f:N_interval ], U_light_01[i_0:i_f:N_interval , args_lyap_modes], cmap="turbo", shading='auto')
    cbar = fig.colorbar(cax_03)
    cbar.ax.tick_params(labelsize=13)
    cbar.set_label('$u(t)$', rotation=0, size=20, labelpad=-50, y=1.1)
    ax03.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
    ax03.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=True)
    ax03.set_xlabel("$i$", fontsize=20)
    ax03.set_ylabel("$t$", fontsize=20)

    cax_04 = ax04.pcolormesh(x_grid, t_light[i_0:i_f:N_interval ], module[i_0:i_f:N_interval, args_lyap_modes], cmap="turbo", shading='auto')
    cbar = fig.colorbar(cax_04)
    cbar.ax.tick_params(labelsize=13)
    cbar.set_label('$r(t)$', rotation=0, size=20, labelpad=-35, y=1.1)
    ax04.tick_params(axis="y", direction="in", labelsize=15, left=False, right=True, labelleft=False, labelright=False)
    ax04.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=True)
    ax04.set_xlabel("$i$", fontsize=20)

    plt.tight_layout()
    plt.savefig(file + '/spatiotemporals.png', dpi=200)
    plt.close()