from back_process import *


if __name__ == "__main__":
    principal_sim = "C:/mnustes_science/simulation_data/FD/PDNLS_extended_PT/extras"
    principal_exp = "C:/mnustes_science/experimental_data"
    dir_05 = principal_sim + "/dimensional/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.0180/sigma=6.000/gamma=0.2800/dist=20.000"
    dir_06 = principal_sim + "/dimensional/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.0180/sigma=6.000/gamma=0.2200/dist=20.000"
    dir_07 = principal_exp + "/PT_05/f=14.20_a=12.00"
    dir_08 = principal_exp + "/PT_05/f=14.20_a=14.00"

    beta = 0.004811649356064012

    Z_r_05 = np.loadtxt(dir_05 + '/field_real.txt', delimiter=',')
    Z_i_05 = np.loadtxt(dir_05 + '/field_img.txt', delimiter=',')
    X_05 = np.loadtxt(dir_05 + '/X.txt', delimiter=',')
    T_05 = np.loadtxt(dir_05 + '/T.txt', delimiter=',')
    Z_complex_05 = Z_r_05 + 1j * Z_i_05
    arg_05 = np.angle(Z_complex_05)
    Z_modulo_05 = np.absolute(Z_complex_05)
    Z_modulo_05 = filtro_superficie(Z_modulo_05, 2, "XY") / np.sqrt(beta)

    ###

    Z_r_06 = np.loadtxt(dir_06 + '/field_real.txt', delimiter=',')
    Z_i_06 = np.loadtxt(dir_06 + '/field_img.txt', delimiter=',')
    X_06 = np.loadtxt(dir_06 + '/X.txt', delimiter=',')
    T_06 = np.loadtxt(dir_06 + '/T.txt', delimiter=',')
    Z_complex_06 = Z_r_06 + 1j * Z_i_06
    arg_06 = np.angle(Z_complex_06)
    Z_modulo_06 = np.absolute(Z_complex_06)
    Z_modulo_06 = filtro_superficie(Z_modulo_06, 2, "XY") / np.sqrt(beta)

    ###

    g = 9790
    l_y = 16
    d = 20
    k_y = np.pi / l_y
    k = k_y
    tau = np.tanh(k * d)
    w_1 = np.sqrt(g * k * tau)

    Z_07 = np.loadtxt(dir_07 + '/Z_strobo.txt', delimiter=',')
    X_07 = np.loadtxt(dir_07 + '/X_mm.txt', delimiter=',')
    T_07 = np.loadtxt(dir_07 + '/T_strobo.txt', delimiter=',')
    injection_07 = np.loadtxt(dir_07 + '/IL_mm.txt', delimiter=',')
    Z_07 = filtro_superficie(Z_07, 5, "YX")

    x07 = injection_07[0] + injection_07[1] - 16 * 2.5
    X_07 = X_07 - x07
    T_07 = T_07 * (w_1 / (2 * np.pi))

    Z_08 = np.loadtxt(dir_08 + '/Z_strobo.txt', delimiter=',')
    X_08 = np.loadtxt(dir_08 + '/X_mm.txt', delimiter=',')
    T_08 = np.loadtxt(dir_08 + '/T_strobo.txt', delimiter=',')
    injection_08 = np.loadtxt(dir_08 + '/IL_mm.txt', delimiter=',')
    Z_08 = filtro_superficie(Z_08, 5, "YX")

    x08 = injection_08[0] + injection_08[1] - 16 * 3 + 2
    X_08 = X_08 - x08
    T_08 = T_08 * (w_1 / (2 * np.pi)) * 2

    xi_exp, xf_exp = -150, 150
    xi_sim, xf_sim = -75, 75
    ti_exp_01, tf_exp_01 = 0, 280
    ti_exp_02, tf_exp_02 = 0, 280
    ti_sim, tf_sim = 0, 280

    data = [[X_08, T_08 - (T_08[-1] - 400), Z_08], [X_07, T_07 - (T_07[-1] - 400), Z_07], [X_05, T_05 - (T_05[-1] - 400), np.flip(Z_modulo_05, axis=1)], [X_06, T_06 - (T_06[-1] - 400), Z_modulo_06]]
    metadata = [[-150, 150], [-150, 150], [-100, 100], [-100, 100]]

    fig, ((ax6, ax5), (ax8, ax7)) = plt.subplots(nrows=2, ncols=2)

    pc_05 = ax5.pcolor(data[0][0], data[0][1], data[0][2], cmap=parula_map)
    ax5.set_ylim(ti_exp_01, tf_exp_01)
    ax5.set_xlim(-150, 150)
    ax5.set_yticklabels([])
    ax5.tick_params(axis="x", direction="in", labeltop=True, labelbottom=False, top=True, bottom=True)
    ax5.tick_params(axis="y", direction="in", left=True, right=True)
    ax5.grid(alpha=0.2, color="k")
    ax5.fill_between([xi_exp, xi_exp + 0.225 * (xf_exp - xi_exp)], [0.705 * (tf_exp_01 - ti_exp_01), tf_exp_01], tf_exp_01, color="white", zorder=12)
    ax5.plot([xi_exp, xi_exp + 0.23 * (xf_exp - xi_exp)], [0.7 * (tf_exp_01 - ti_exp_01), tf_exp_01], linewidth=1, color="k", zorder=12)
    ax5.plot([xi_exp, xi_exp + 0.23 * (xf_exp - xi_exp)], [tf_exp_01, tf_exp_01], linewidth=1, color="k", zorder=12)
    ax5.plot([xi_exp, xi_exp], [0.7 * (tf_exp_01 - ti_exp_01), tf_exp_01], linewidth=1, color="k", zorder=12)
    ax5.text(xi_exp + 0.005 * (xf_exp - xi_exp), ti_exp_01 + 0.795 * (tf_exp_01 - ti_exp_01), "$\\textbf{EXP}$", fontsize=13, rotation=45, zorder=13)
    ax5.text(xi_exp + 0.84 * (xf_exp - xi_exp), ti_exp_01 + 0.87 * (tf_exp_01 - ti_exp_01), "$\\textbf{(c)}$", fontsize=18, rotation=0, zorder=13, color="white")

    pc_06 = ax6.pcolor(data[1][0], data[1][1], data[1][2], cmap=parula_map)
    ax6.set_ylim(ti_exp_02, tf_exp_02)
    ax6.set_xlim(-150, 150)
    ax6.tick_params(axis="x", direction="in", labeltop=True, labelbottom=False, top=True, bottom=True)
    ax6.tick_params(axis="y", direction="in", left=True, right=True)
    ax6.grid(alpha=0.2, color="k")
    ax6.set_yticklabels([])
    ax6.fill_between([xi_exp, xi_exp + 0.225 * (xf_exp - xi_exp)], [0.705 * (tf_exp_02 - ti_exp_02), tf_exp_02], tf_exp_02, color="white", zorder=12)
    ax6.plot([xi_exp, xi_exp + 0.23 * (xf_exp - xi_exp)], [0.7 * (tf_exp_02 - ti_exp_02), tf_exp_02], linewidth=1, color="k", zorder=12)
    ax6.plot([xi_exp, xi_exp + 0.23 * (xf_exp - xi_exp)], [tf_exp_02, tf_exp_02], linewidth=1, color="k", zorder=12)
    ax6.plot([xi_exp, xi_exp], [0.7 * (tf_exp_02 - ti_exp_02), tf_exp_02], linewidth=1, color="k", zorder=12)
    ax6.text(xi_exp + 0.005 * (xf_exp - xi_exp), ti_exp_02 + 0.795 * (tf_exp_02 - ti_exp_02), "$\\textbf{EXP}$", fontsize=13, rotation=45, zorder=13)
    ax6.text(xi_exp + 0.84 * (xf_exp - xi_exp), ti_exp_02 + 0.87 * (tf_exp_02 - ti_exp_02), "$\\textbf{(b)}$", fontsize=18, rotation=0, zorder=13, color="white")

    pc_07 = ax7.pcolor(data[2][0], data[2][1], data[2][2], cmap=parula_map)
    ax7.set_ylim(ti_sim, tf_sim)
    ax7.set_xlim(xi_sim, xf_sim)
    ax7.grid(alpha=0.2, color="k")
    ax7.tick_params(axis="x", direction="in", top=True, bottom=True)
    ax7.tick_params(axis="y", direction="in", left=True, right=True)
    ax7.set_yticklabels([])
    ax7.fill_between([xi_sim, xi_sim + 0.225 * (xf_sim - xi_sim)], [0.705 * (tf_sim - ti_sim), tf_sim], tf_sim, color="white", zorder=12)
    ax7.plot([xi_sim, xi_sim + 0.23 * (xf_sim - xi_sim)], [0.7 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax7.plot([xi_sim, xi_sim + 0.23 * (xf_sim - xi_sim)], [tf_sim, tf_sim], linewidth=1, color="k", zorder=12)
    ax7.plot([xi_sim, xi_sim], [0.7 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax7.text(xi_sim + 0.005 * (xf_sim - xi_sim), ti_sim + 0.795 * (tf_sim - ti_sim), "$\\textbf{NUM}$", fontsize=13, rotation=45, zorder=13)
    ax7.text(xi_sim + 0.84 * (xf_sim - xi_sim), ti_sim + 0.87 * (tf_sim - ti_sim), "$\\textbf{(f)}$", fontsize=18, rotation=0, zorder=13, color="white")

    pc_08 = ax8.pcolor(data[3][0], data[3][1], data[3][2], cmap=parula_map)
    ax8.set_ylim(ti_sim, tf_sim)
    ax8.set_xlim(xi_sim, xf_sim)
    ax8.grid(alpha=0.1, color="k")
    ax8.tick_params(axis="x", direction="in", top=True, bottom=True)
    ax8.tick_params(axis="y", direction="in", left=True, right=True, zorder=10)
    ax8.set_yticklabels([])
    ax8.fill_between([xi_sim, xi_sim + 0.225 * (xf_sim - xi_sim)], [0.705 * (tf_sim - ti_sim), tf_sim], tf_sim, color="white", zorder=12)
    ax8.plot([xi_sim, xi_sim + 0.23 * (xf_sim - xi_sim)], [0.7 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax8.plot([xi_sim, xi_sim + 0.23 * (xf_sim - xi_sim)], [tf_sim, tf_sim], linewidth=1, color="k", zorder=12)
    ax8.plot([xi_sim, xi_sim], [0.7 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax8.text(xi_sim + 0.005 * (xf_sim - xi_sim), ti_sim + 0.795 * (tf_sim - ti_sim), "$\\textbf{NUM}$", fontsize=13, rotation=45, zorder=13)
    ax8.text(xi_sim + 0.84 * (xf_sim - xi_sim), ti_sim + 0.87 * (tf_sim - ti_sim), "$\\textbf{(e)}$", fontsize=18, rotation=0, zorder=13, color="white")


    plt.subplots_adjust(wspace=0.2, hspace=0.07)
    cax_01 = fig.add_axes([0.49, 0.51, 0.015, 0.37])
    fig.colorbar(pc_05, cax=cax_01)
    cax_02 = fig.add_axes([0.913, 0.51, 0.015, 0.37])
    fig.colorbar(pc_06, cax=cax_02)
    cax_03 = fig.add_axes([0.49, 0.11, 0.015, 0.37])
    fig.colorbar(pc_07, cax=cax_03)
    cax_04 = fig.add_axes([0.913, 0.11, 0.015, 0.37])
    fig.colorbar(pc_08, cax=cax_04)
    savedirectory = "C:/mnustes_science/simulation_data/FD/PDNLS_extended_PT/extras/dimensional/analysis/figures"
    plt.savefig('Fig04.png', dpi=300)
    plt.close()