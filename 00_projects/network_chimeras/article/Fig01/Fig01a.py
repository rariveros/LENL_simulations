from back_process import *

if __name__ == '__main__':

    dir_01 = "D:/mnustes_science/simulation_data/FD/network_chimeras/FIG01/erdos_renyi/k=0.0150/mean_degree=18.00"
    dir_02 = "D:/mnustes_science/simulation_data/FD/network_chimeras/FIG01/erdos_renyi/k=0.0150/mean_degree=22.00"
    dir_03 = "D:/mnustes_science/simulation_data/FD/network_chimeras/FIG01/erdos_renyi/k=0.0150/mean_degree=26.00"

    U_light_01 = np.loadtxt(dir_01 + '/U_lyap.txt', delimiter=',')
    V_light_01 = np.loadtxt(dir_01 + '/V_lyap.txt', delimiter=',')
    arg_lyap_pos_01 = np.loadtxt(dir_01 + '/arg_lyap.txt', delimiter=',')
    x_01 = np.loadtxt(dir_01 + '/x_lyap.txt', delimiter=',')
    lyap_mean_01 = np.loadtxt(dir_01 + '/lyap_spectrum.txt', delimiter=',')
    Q_pos_01 = np.loadtxt(dir_01 + '/Q_pos.txt', delimiter=',')

    U_light_02 = np.loadtxt(dir_02 + '/U_lyap.txt', delimiter=',')
    V_light_02 = np.loadtxt(dir_02 + '/V_lyap.txt', delimiter=',')
    arg_lyap_02 = np.loadtxt(dir_02 + '/arg_lyap.txt', delimiter=',')
    x_02 = np.loadtxt(dir_02 + '/x_lyap.txt', delimiter=',')
    lyap_mean_02 = np.loadtxt(dir_02 + '/lyap_spectrum.txt', delimiter=',')
    Q_pos_02 = np.loadtxt(dir_02 + '/Q_pos.txt', delimiter=',')

    U_light_03 = np.loadtxt(dir_03 + '/U_lyap.txt', delimiter=',')
    V_light_03 = np.loadtxt(dir_03 + '/V_lyap.txt', delimiter=',')
    arg_lyap_03 = np.loadtxt(dir_03 + '/arg_lyap.txt', delimiter=',')
    x_03 = np.loadtxt(dir_03 + '/x_lyap.txt', delimiter=',')
    lyap_mean_03 = np.loadtxt(dir_03 + '/lyap_spectrum.txt', delimiter=',')
    Q_pos_03 = np.loadtxt(dir_03 + '/Q_pos.txt', delimiter=',')

    t_light = np.loadtxt(dir_01 + '/T_lyap.txt', delimiter=',')
    x_grid = np.loadtxt(dir_01 + '/X.txt', delimiter=',')

    module_01 = np.sqrt(U_light_01 ** 2 + V_light_01 ** 2)
    module_02 = np.sqrt(U_light_02 ** 2 + V_light_02 ** 2)
    module_03 = np.sqrt(U_light_03 ** 2 + V_light_03 ** 2)

    dt = t_light[1] - t_light[0]
    Nt = len(t_light)
    Nx = len(x_grid)

    t_init = 4980
    t_final = 5000
    i_0 = np.argmin(np.abs(t_light - t_init))
    i_f = np.argmin(np.abs(t_light - t_final))
    power_threshold = 1.5

    file = "D:/mnustes_science/simulation_data/FD/network_chimeras/FIG01/erdos_renyi/k=0.0150"

    arg_lyap_pos_01 = lyap_mean_01 >= 0.0
    args_lyap_modes_01 = np.argsort(Q_pos_01)

    arg_lyap_pos_02 = lyap_mean_02 >= 0.0
    args_lyap_modes_02 = np.argsort(Q_pos_02)

    arg_lyap_pos_03 = lyap_mean_03 >= 0.0
    args_lyap_modes_03 = np.argsort(Q_pos_03)

    #### SPATIOTEMPORAL DIAGRAMS ####
    N_interval = 2
    fig, (ax01, ax02, ax03) = plt.subplots(3, 1, figsize=(2.5, 4))
    cax_01 = ax01.pcolormesh(x_grid, t_light[i_0:i_f:N_interval] + 1, module_01[i_0:i_f:N_interval, :], cmap="turbo", shading='auto')
    cbar = fig.colorbar(cax_01, pad=0.03, aspect=10)
    cbar.ax.tick_params(labelsize=13)
    cbar.set_label('$r(t)$', rotation=0, size=18, labelpad=-10, y=1.3)
    ax01.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
    ax01.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=False, labeltop=False, labelbottom=False)

    cax_02 = ax02.pcolormesh(x_grid, t_light[i_0:i_f:N_interval ] + 1, module_02[i_0:i_f:N_interval, :], cmap="turbo", shading='auto')
    cbar = fig.colorbar(cax_02, pad=0.03, aspect=10)
    cbar.ax.tick_params(labelsize=13)
    #cbar.set_label('$r(t)$', rotation=0, size=20, labelpad=-50, y=1.1)
    ax02.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
    ax02.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=False, labeltop=False, labelbottom=False)
    ax02.set_ylabel("$\\textrm{Time}$", fontsize=18)

    cax_03 = ax03.pcolormesh(x_grid, t_light[i_0:i_f:N_interval ] + 1, module_03[i_0:i_f:N_interval, :], cmap="turbo", shading='auto')
    cbar = fig.colorbar(cax_03, pad=0.03, aspect=10)
    cbar.ax.tick_params(labelsize=13)
    #cbar.set_label('$r(t)$', rotation=0, size=20, labelpad=-50, y=1.1)
    ax03.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
    ax03.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=True)
    ax03.set_xlabel("$\\textrm{Index}$", fontsize=15)

    plt.subplots_adjust(
        left=0.12,
        right=0.88,
        top=0.95,
        bottom=0.07,
        hspace=-0.03,
        wspace=0.0# controls vertical spacing
    )
    plt.tight_layout()
    plt.savefig(file + '/spatiotemporals.png', dpi=200)
    plt.close()

    #### LYAPUNOV SPECTRUM ####
    N_interval = 2
    fig, ax = plt.subplots(1, 1, figsize=(3.5, 4))
    ax.plot(x_01, lyap_mean_01, c="k", linewidth=2, zorder=0)
    ax.scatter(x_01, lyap_mean_01, color="r", edgecolors="k", linewidth=0.8, zorder=2, s=12)
    ax.plot(x_02, lyap_mean_02, c="k", linewidth=2, zorder=0)
    ax.scatter(x_02, lyap_mean_02, color="g", edgecolors="k", linewidth=0.8, zorder=2, s=12)
    ax.plot(x_03, lyap_mean_03, c="k", linewidth=2, zorder=0)
    ax.scatter(x_03, lyap_mean_03, color="b", edgecolors="k", linewidth=0.8, zorder=2, s=12)

    plt.hlines(0, 1, 0, color="k", linestyle="-", zorder=0)
    ax.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
    ax.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=True)
    ax.set_xlabel("$\\textrm{Normalized index}$", fontsize=15)
    ax.set_ylabel("$\\textrm{Lyapunov spectrum}$", fontsize=15)
    ax.set_xlim(0, 0.5)
    ax.set_ylim(-0.05, 0.12)

    plt.subplots_adjust(
        left=0.12,
        right=0.88,
        top=0.95,
        bottom=0.1,
        hspace=-0.03,
        wspace=0.0# controls vertical spacing
    )

    plt.tight_layout()
    plt.savefig(file + '/lyap_spectrums.png', dpi=200)
    plt.close()