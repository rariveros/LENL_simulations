from back_process import *


if __name__ == "__main__":
    principal_sim = "C:/mnustes_science/simulation_data/FD/PDNLS_extended_PT/extras"
    principal_exp = "C:/mnustes_science/experimental_data"
    dir_01 = principal_sim + "/dist_03/alpha=13.0480/beta=1.000/mu=0.0950/nu=0.0180/sigma=12.6700/gamma=0.1690/dist=40.000"
    #"/dimensional_dist_02/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.0180/sigma=6.000/gamma=0.1850/dist=31.000"
    dir_02 = principal_sim + "/dist_03/alpha=13.0480/beta=1.000/mu=0.0950/nu=0.0180/sigma=12.6700/gamma=0.1690/dist=31.000"
    #"/dimensional/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.0180/sigma=6.000/gamma=0.1850/dist=20.000"
    dir_03 = principal_exp + "/PT_04/f=14.20_a=10.00"
    dir_04 = principal_exp + "/PT_06/f=14.20_a=10.50_d2"

    beta = 0.004811649356064012

    Z_r_01 = np.loadtxt(dir_01 + '/field_real.txt', delimiter=',')
    Z_i_01 = np.loadtxt(dir_01 + '/field_img.txt', delimiter=',')
    X_01 = np.loadtxt(dir_01 + '/X.txt', delimiter=',')
    T_01 = np.loadtxt(dir_01 + '/T.txt', delimiter=',')
    Z_complex_01 = Z_r_01 + 1j * Z_i_01
    arg_01 = np.angle(Z_complex_01)
    Z_modulo_01 = np.absolute(Z_complex_01)
    Z_modulo_01 = filtro_superficie(Z_modulo_01, 5, "XY") / np.sqrt(beta)

    ###

    Z_r_02 = np.loadtxt(dir_02 + '/field_real.txt', delimiter=',')
    Z_i_02 = np.loadtxt(dir_02 + '/field_img.txt', delimiter=',')
    X_02 = np.loadtxt(dir_02 + '/X.txt', delimiter=',')
    T_02 = np.loadtxt(dir_02 + '/T.txt', delimiter=',')
    Z_complex_02 = Z_r_02 + 1j * Z_i_02
    arg_02 = np.angle(Z_complex_02)
    Z_modulo_02 = np.absolute(Z_complex_02)
    Z_modulo_02 = filtro_superficie(Z_modulo_02, 5, "XY") / np.sqrt(beta)

    ###

    g = 9790
    l_y = 16
    d = 20
    k_y = np.pi / l_y
    k = k_y
    tau = np.tanh(k * d)
    w_1 = np.sqrt(g * k * tau)

    Z_03 = np.loadtxt(dir_03 + '/Z_strobo.txt', delimiter=',')
    X_03 = np.loadtxt(dir_03 + '/X_mm.txt', delimiter=',')
    T_03 = np.loadtxt(dir_03 + '/T_strobo.txt', delimiter=',')
    injection_03 = np.loadtxt(dir_03 + '/IL_mm.txt', delimiter=',')
    Z_03 = filtro_superficie(Z_03, 5, "YX")

    x03 = -9
    X_03 = X_03 - x03
    T_03 = T_03 * (w_1 / (2 * np.pi))

    Z_04 = np.loadtxt(dir_04 + '/Z_strobo.txt', delimiter=',')
    X_04 = np.loadtxt(dir_04 + '/X_mm.txt', delimiter=',')
    T_04 = np.loadtxt(dir_04 + '/T_strobo.txt', delimiter=',')
    injection_04 = np.loadtxt(dir_04 + '/IL_mm.txt', delimiter=',')
    Z_04 = filtro_superficie(Z_04, 5, "YX")
    x04 = injection_04[0] + injection_04[1]
    X_04 = X_04 - x04
    T_04 = T_04 * (w_1 / (2 * np.pi))

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
    Z_modulo_05 = filtro_superficie(Z_modulo_05, 5, "XY") / np.sqrt(beta)

    ###

    Z_r_06 = np.loadtxt(dir_06 + '/field_real.txt', delimiter=',')
    Z_i_06 = np.loadtxt(dir_06 + '/field_img.txt', delimiter=',')
    X_06 = np.loadtxt(dir_06 + '/X.txt', delimiter=',')
    T_06 = np.loadtxt(dir_06 + '/T.txt', delimiter=',')
    Z_complex_06 = Z_r_06 + 1j * Z_i_06
    arg_06 = np.angle(Z_complex_06)
    Z_modulo_06 = np.absolute(Z_complex_06)
    Z_modulo_06 = filtro_superficie(Z_modulo_06, 5, "XY") / np.sqrt(beta)

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

    x08 = injection_08[0] + injection_08[1] - 16 * 3
    X_08 = X_08 - x08
    T_08 = T_08 * (w_1 / (2 * np.pi))

    xi_exp, xf_exp = -150, 150
    xi_sim, xf_sim = -150, 150#-150/np.sqrt(2), 150/np.sqrt(2)
    ti_exp_01, tf_exp_01 = 0, 400
    ti_exp_02, tf_exp_02 = 0, 400
    ti_sim, tf_sim = 0, 400

    data = [[X_04, T_04 - (T_04[-1] - 400), Z_04], [X_03, T_03 - (T_03[-1] - 400), Z_03], [X_01, T_01 - (T_01[-1] - 400), np.flip(Z_modulo_01, axis=1)], [X_02, T_02 - (T_02[-1] - 400), Z_modulo_02]]
    metadata = [[-150, 150], [-150, 150], [-100, 100], [-100, 100]]

    fig, ((ax1, ax2, ax5, ax6), (ax3, ax4, ax7, ax8)) = plt.subplots(nrows=2, ncols=4, figsize=(16, 4))

    pc_01 = ax1.pcolor(data[0][0], data[0][1], data[0][2], cmap=parula_map)
    ax1.set_ylim(0, 400)
    ax1.set_xlim(-150, 150)
    ax1.tick_params(axis="x", direction="in", labeltop=True, labelbottom=False, top=True, bottom=True)
    ax1.tick_params(axis="y", direction="in", left=True, right=True)
    ax1.grid(alpha=0.2, color="k")
    ax1.fill_between([xi_exp, xi_exp + 0.225 * (xf_exp - xi_exp)], [0.705 * (tf_exp_01 - ti_exp_01), tf_exp_01], tf_exp_01, color="white", zorder=12)
    ax1.plot([xi_exp, xi_exp + 0.23 * (xf_exp - xi_exp)], [0.7 * (tf_exp_01 - ti_exp_01), tf_exp_01], linewidth=1, color="k", zorder=12)
    ax1.plot([xi_exp, xi_exp + 0.23 * (xf_exp - xi_exp)], [tf_exp_01, tf_exp_01], linewidth=1, color="k", zorder=12)
    ax1.plot([xi_exp, xi_exp], [0.7 * (tf_exp_01 - ti_exp_01), tf_exp_01], linewidth=1, color="k", zorder=12)
    ax1.text(xi_exp + 0.005 * (xf_exp - xi_exp), ti_exp_01 + 0.795 * (tf_exp_01 - ti_exp_01), "$\\textbf{EXP}$", fontsize=13, rotation=45, zorder=13)

    pc_02 = ax2.pcolor(data[1][0], data[1][1], data[1][2], cmap=parula_map)
    ax2.set_ylim(0, 400)
    ax2.set_yticklabels([])
    ax2.set_xlim(-150, 150)
    ax2.tick_params(axis="x", direction="in", labeltop=True, labelbottom=False, top=True, bottom=True)
    ax2.tick_params(axis="y", direction="in", left=True, right=True)
    ax2.grid(alpha=0.2, color="k")
    ax2.fill_between([xi_exp, xi_exp + 0.225 * (xf_exp - xi_exp)], [0.705 * (tf_exp_02 - ti_exp_02), tf_exp_02], tf_exp_02, color="white", zorder=12)
    ax2.plot([xi_exp, xi_exp + 0.23 * (xf_exp - xi_exp)], [0.7 * (tf_exp_02 - ti_exp_02), tf_exp_02], linewidth=1, color="k", zorder=12)
    ax2.plot([xi_exp, xi_exp + 0.23 * (xf_exp - xi_exp)], [tf_exp_02, tf_exp_02], linewidth=1, color="k", zorder=12)
    ax2.plot([xi_exp, xi_exp], [0.7 * (tf_exp_02 - ti_exp_02), tf_exp_02], linewidth=1, color="k", zorder=12)
    ax2.text(xi_exp + 0.005 * (xf_exp - xi_exp), ti_exp_02 + 0.795 * (tf_exp_02 - ti_exp_02), "$\\textbf{EXP}$", fontsize=13, rotation=45, zorder=13)

    pc_03 = ax3.pcolor(data[2][0], data[2][1], data[2][2], cmap=parula_map)
    ax3.set_ylim(0, 400)
    ax3.set_xlim(xi_sim, xf_sim)
    ax3.grid(alpha=0.2, color="k")
    ax3.tick_params(axis="x", direction="in", top=True, bottom=True)
    ax3.tick_params(axis="y", direction="in", left=True, right=True)
    ax3.fill_between([xi_sim, xi_sim + 0.225 * (xf_sim - xi_sim)], [0.705 * (tf_sim - ti_sim), tf_sim], tf_sim, color="white", zorder=12)
    ax3.plot([xi_sim, xi_sim + 0.23 * (xf_sim - xi_sim)], [0.7 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax3.plot([xi_sim, xi_sim + 0.23 * (xf_sim - xi_sim)], [tf_sim, tf_sim], linewidth=1, color="k", zorder=12)
    ax3.plot([xi_sim, xi_sim], [0.7 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax3.text(xi_sim + 0.005 * (xf_sim - xi_sim), ti_sim + 0.795 * (tf_sim - ti_sim), "$\\textbf{NUM}$", fontsize=13, rotation=45, zorder=13)

    pc_04 = ax4.pcolor(data[3][0], data[3][1], data[3][2], cmap=parula_map)
    ax4.set_ylim(0, 400)
    ax4.set_yticklabels([])
    ax4.set_xlim(xi_sim, xf_sim)
    ax4.grid(alpha=0.1, color="k")
    ax4.tick_params(axis="x", direction="in", top=True, bottom=True)
    ax4.tick_params(axis="y", direction="in", left=True, right=True, zorder=10)
    ax4.fill_between([xi_sim, xi_sim + 0.225 * (xf_sim - xi_sim)], [0.705 * (tf_sim - ti_sim), tf_sim], tf_sim, color="white", zorder=12)
    ax4.plot([xi_sim, xi_sim + 0.23 * (xf_sim - xi_sim)], [0.7 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax4.plot([xi_sim, xi_sim + 0.23 * (xf_sim - xi_sim)], [tf_sim, tf_sim], linewidth=1, color="k", zorder=12)
    ax4.plot([xi_sim, xi_sim], [0.7 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax4.text(xi_sim + 0.005 * (xf_sim - xi_sim), ti_sim + 0.795 * (tf_sim - ti_sim), "$\\textbf{NUM}$", fontsize=13, rotation=45, zorder=13)

    data = [[X_08, T_08 - (T_08[-1] - 400), Z_08], [X_07, T_07 - (T_07[-1] - 400), Z_07],
            [X_05, T_05 - (T_05[-1] - 400), np.flip(Z_modulo_05, axis=1)], [X_06, T_06 - (T_06[-1] - 400), Z_modulo_06]]
    metadata = [[-150, 150], [-150, 150], [-100, 100], [-100, 100]]


    ax5.pcolor(data[0][0], data[0][1], data[0][2], cmap=parula_map)
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

    ax6.pcolor(data[1][0], data[1][1], data[1][2], cmap=parula_map)
    ax6.set_ylim(ti_exp_02, tf_exp_02)

    ax6.set_xlim(-150, 150)
    ax6.tick_params(axis="x", direction="in", labeltop=True, labelbottom=False, top=True, bottom=True)
    ax6.tick_params(axis="y", direction="in", left=True, right=True)
    ax6.grid(alpha=0.2, color="k")
    ax6.fill_between([xi_exp, xi_exp + 0.225 * (xf_exp - xi_exp)], [0.705 * (tf_exp_02 - ti_exp_02), tf_exp_02], tf_exp_02, color="white", zorder=12)
    ax6.plot([xi_exp, xi_exp + 0.23 * (xf_exp - xi_exp)], [0.7 * (tf_exp_02 - ti_exp_02), tf_exp_02], linewidth=1, color="k", zorder=12)
    ax6.plot([xi_exp, xi_exp + 0.23 * (xf_exp - xi_exp)], [tf_exp_02, tf_exp_02], linewidth=1, color="k", zorder=12)
    ax6.plot([xi_exp, xi_exp], [0.7 * (tf_exp_02 - ti_exp_02), tf_exp_02], linewidth=1, color="k", zorder=12)
    ax6.text(xi_exp + 0.005 * (xf_exp - xi_exp), ti_exp_02 + 0.795 * (tf_exp_02 - ti_exp_02), "$\\textbf{EXP}$", fontsize=13, rotation=45, zorder=13)

    ax7.pcolor(data[2][0], data[2][1], data[2][2], cmap=parula_map)
    ax7.set_ylim(0, 400)
    ax7.set_xlim(-90, 90)
    ax7.grid(alpha=0.2, color="k")
    ax7.tick_params(axis="x", direction="in", top=True, bottom=True)
    ax7.tick_params(axis="y", direction="in", left=True, right=True)
    ax7.set_yticklabels([])
    ax7.fill_between([xi_sim, xi_sim + 0.225 * (xf_sim - xi_sim)], [0.705 * (tf_sim - ti_sim), tf_sim], tf_sim, color="white", zorder=12)
    ax7.plot([xi_sim, xi_sim + 0.23 * (xf_sim - xi_sim)], [0.7 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax7.plot([xi_sim, xi_sim + 0.23 * (xf_sim - xi_sim)], [tf_sim, tf_sim], linewidth=1, color="k", zorder=12)
    ax7.plot([xi_sim, xi_sim], [0.7 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax7.text(xi_sim + 0.005 * (xf_sim - xi_sim), ti_sim + 0.795 * (tf_sim - ti_sim), "$\\textbf{NUM}$", fontsize=13, rotation=45, zorder=13)

    pc = ax8.pcolor(data[3][0], data[3][1], data[3][2], cmap=parula_map)
    ax8.set_ylim(0, 400)
    ax8.set_xlim(-90, 90)
    ax8.grid(alpha=0.1, color="k")
    ax8.tick_params(axis="x", direction="in", top=True, bottom=True)
    ax8.tick_params(axis="y", direction="in", left=True, right=True, zorder=10)
    ax8.fill_between([xi_sim, xi_sim + 0.225 * (xf_sim - xi_sim)], [0.705 * (tf_sim - ti_sim), tf_sim], tf_sim, color="white", zorder=12)
    ax8.plot([xi_sim, xi_sim + 0.23 * (xf_sim - xi_sim)], [0.7 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax8.plot([xi_sim, xi_sim + 0.23 * (xf_sim - xi_sim)], [tf_sim, tf_sim], linewidth=1, color="k", zorder=12)
    ax8.plot([xi_sim, xi_sim], [0.7 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax8.text(xi_sim + 0.005 * (xf_sim - xi_sim), ti_sim + 0.795 * (tf_sim - ti_sim), "$\\textbf{NUM}$", fontsize=13, rotation=45, zorder=13)


    plt.subplots_adjust(wspace=0.2, hspace=0.07)
    cax_01 = fig.add_axes([0.49, 0.51, 0.015, 0.37])
    fig.colorbar(pc_01, cax=cax_01)
    cax_02 = fig.add_axes([0.913, 0.51, 0.015, 0.37])
    fig.colorbar(pc_02, cax=cax_02)
    cax_03 = fig.add_axes([0.49, 0.11, 0.015, 0.37])
    fig.colorbar(pc_03, cax=cax_03)
    cax_04 = fig.add_axes([0.913, 0.11, 0.015, 0.37])
    fig.colorbar(pc_04, cax=cax_04)
    savedirectory = "C:/mnustes_science/simulation_data/FD/PDNLS_extended_PT/extras/dimensional/analysis/figures"
    plt.savefig(savedirectory + '/Fig02_poster.png', dpi=300)
    plt.close()