from back_process import *


if __name__ == "__main__":
    principal_sim = "E:/mnustes_science/simulation_data/FD/PDNLS_extended_PT/extras"
    principal_exp = "E:/mnustes_science/experimental_data"
    dir_01 = principal_sim + "/dimensional_dist_02/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.0180/sigma=6.000/gamma=0.1850/dist=31.000"#/dist_03/alpha=13.0480/beta=1.000/mu=0.0950/nu=0.0180/sigma=12.6700/gamma=0.1690/dist=40.000"#
    dir_02 = principal_sim + "/dimensional_dist_02/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.0180/sigma=6.000/gamma=0.1850/dist=21.0000"#"/dist_03/alpha=13.0480/beta=1.000/mu=0.0950/nu=0.0180/sigma=12.6700/gamma=0.1690/dist=31.000"#
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

    xi_exp, xf_exp = -150, 150
    xi_sim, xf_sim = -75, 75
    ti_exp_01, tf_exp_01 = 0, 400
    ti_exp_02, tf_exp_02 = 0, 400
    ti_sim, tf_sim = 0, 400

    data = [[X_04, T_04 - (T_04[-1] - 400), Z_04], [X_03, T_03 - (T_03[-1] - 400), Z_03], [X_01, T_01 - (T_01[-1] - 400), np.flip(Z_modulo_01, axis=1)], [X_02, T_02 - (T_02[-1] - 400), Z_modulo_02]]
    metadata = [[-150, 150], [-150, 150], [-100, 100], [-100, 100]]

    fig, ((ax2, ax1), (ax4, ax3)) = plt.subplots(nrows=2, ncols=2)

    pc_01 = ax1.pcolor(data[0][0], data[0][1], data[0][2], cmap=parula_map)
    ax1.set_ylim(0, 400)
    ax1.set_xlim(-150, 150)

    ax1.set_yticklabels([])
    ax1.tick_params(axis="x", direction="in", labeltop=True, labelbottom=False, top=True, bottom=True, labelsize=12)
    ax1.tick_params(axis="y", direction="in", left=True, right=True, labelsize=12)
    ax1.grid(alpha=0.2, color="k")

    ax1.fill_between([xi_exp, xi_exp + 0.225 * (xf_exp - xi_exp)], [0.705 * (tf_exp_01 - ti_exp_01), tf_exp_01], tf_exp_01, color="white", zorder=12)
    ax1.plot([xi_exp, xi_exp + 0.23 * (xf_exp - xi_exp)], [0.7 * (tf_exp_01 - ti_exp_01), tf_exp_01], linewidth=1, color="k", zorder=12)
    ax1.plot([xi_exp, xi_exp + 0.23 * (xf_exp - xi_exp)], [tf_exp_01, tf_exp_01], linewidth=1, color="k", zorder=12)
    ax1.plot([xi_exp, xi_exp], [0.7 * (tf_exp_01 - ti_exp_01), tf_exp_01], linewidth=1, color="k", zorder=12)
    ax1.text(xi_exp + 0.005 * (xf_exp - xi_exp), ti_exp_01 + 0.795 * (tf_exp_01 - ti_exp_01), "$\\textbf{EXP}$", fontsize=13, rotation=45, zorder=13)
    ax1.text(xi_exp + 0.84 * (xf_exp - xi_exp), ti_exp_01 + 0.87 * (tf_exp_01 - ti_exp_01), "$\\textbf{(b)}$", fontsize=18, rotation=0, zorder=13, color="white")

    pc_02 = ax2.pcolor(data[1][0], data[1][1], data[1][2], cmap=parula_map)
    ax2.set_ylim(0, 400)
    ax2.set_xlim(-150, 150)

    ax2.tick_params(axis="x", direction="in", labeltop=True, labelbottom=False, top=True, bottom=True, labelsize=12)
    ax2.tick_params(axis="y", direction="in", left=True, right=True, labelsize=12)
    ax2.set_ylabel("$t/T$", fontsize=15)

    ax2.grid(alpha=0.2, color="k")
    ax2.fill_between([xi_exp, xi_exp + 0.225 * (xf_exp - xi_exp)], [0.705 * (tf_exp_02 - ti_exp_02), tf_exp_02], tf_exp_02, color="white", zorder=12)
    ax2.plot([xi_exp, xi_exp + 0.23 * (xf_exp - xi_exp)], [0.7 * (tf_exp_02 - ti_exp_02), tf_exp_02], linewidth=1, color="k", zorder=12)
    ax2.plot([xi_exp, xi_exp + 0.23 * (xf_exp - xi_exp)], [tf_exp_02, tf_exp_02], linewidth=1, color="k", zorder=12)
    ax2.plot([xi_exp, xi_exp], [0.7 * (tf_exp_02 - ti_exp_02), tf_exp_02], linewidth=1, color="k", zorder=12)
    ax2.text(xi_exp + 0.005 * (xf_exp - xi_exp), ti_exp_02 + 0.795 * (tf_exp_02 - ti_exp_02), "$\\textbf{EXP}$", fontsize=13, rotation=45, zorder=13)
    ax2.text(xi_exp + 0.84 * (xf_exp - xi_exp), ti_exp_02 + 0.87 * (tf_exp_02 - ti_exp_02), "$\\textbf{(a)}$", fontsize=18, rotation=0, zorder=13, color="white")

    pc_03 = ax3.pcolor(data[2][0], data[2][1], data[2][2], cmap=parula_map)
    ax3.set_ylim(0, 400)
    ax3.set_xlim(xi_sim, xf_sim)
    ax3.grid(alpha=0.2, color="k")

    ax3.set_yticklabels([])
    ax3.set_xlabel("$x\ \\textrm{(mm)}$", fontsize=15)
    ax3.tick_params(axis="x", direction="in", top=True, bottom=True, labelsize=12)
    ax3.tick_params(axis="y", direction="in", left=True, right=True, labelsize=12)

    ax3.fill_between([xi_sim, xi_sim + 0.225 * (xf_sim - xi_sim)], [0.705 * (tf_sim - ti_sim), tf_sim], tf_sim, color="white", zorder=12)
    ax3.plot([xi_sim, xi_sim + 0.23 * (xf_sim - xi_sim)], [0.7 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax3.plot([xi_sim, xi_sim + 0.23 * (xf_sim - xi_sim)], [tf_sim, tf_sim], linewidth=1, color="k", zorder=12)
    ax3.plot([xi_sim, xi_sim], [0.7 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax3.text(xi_sim + 0.005 * (xf_sim - xi_sim), ti_sim + 0.795 * (tf_sim - ti_sim), "$\\textbf{NUM}$", fontsize=13, rotation=45, zorder=13)
    ax3.text(xi_sim + 0.84 * (xf_sim - xi_sim), ti_sim + 0.87 * (tf_sim - ti_sim), "$\\textbf{(d)}$", fontsize=18, rotation=0, zorder=13, color="white")

    pc_04 = ax4.pcolor(data[3][0], data[3][1] + 0, data[3][2], cmap=parula_map)
    ax4.set_ylim(0, 400)
    ax4.set_xlim(xi_sim, xf_sim)
    ax4.grid(alpha=0.2, color="k")

    ax4.set_xlabel("$x\ \\textrm{(mm)}$", fontsize=15)
    ax4.set_ylabel("$t/T$", fontsize=15)
    ax4.tick_params(axis="x", direction="in", top=True, bottom=True, labelsize=12)
    ax4.tick_params(axis="y", direction="in", left=True, right=True, zorder=10, labelsize=12)

    ax4.fill_between([xi_sim, xi_sim + 0.225 * (xf_sim - xi_sim)], [0.705 * (tf_sim - ti_sim), tf_sim], tf_sim, color="white", zorder=12)
    ax4.plot([xi_sim, xi_sim + 0.23 * (xf_sim - xi_sim)], [0.7 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax4.plot([xi_sim, xi_sim + 0.23 * (xf_sim - xi_sim)], [tf_sim, tf_sim], linewidth=1, color="k", zorder=12)
    ax4.plot([xi_sim, xi_sim], [0.7 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax4.text(xi_sim + 0.005 * (xf_sim - xi_sim), ti_sim + 0.795 * (tf_sim - ti_sim), "$\\textbf{NUM}$", fontsize=13, rotation=45, zorder=13)
    ax4.text(xi_sim + 0.84 * (xf_sim - xi_sim), ti_sim + 0.87 * (tf_sim - ti_sim), "$\\textbf{(c)}$", fontsize=18, rotation=0, zorder=13, color="white")

    plt.subplots_adjust(wspace=0.2, hspace=0.07)
    cax_02 = fig.add_axes([0.49, 0.51, 0.015, 0.37])
    cbar_02 = fig.colorbar(pc_02, cax=cax_02)
    cbar_02.ax.tick_params(labelsize=11)

    cax_01 = fig.add_axes([0.913, 0.51, 0.015, 0.37])
    cbar_01 = fig.colorbar(pc_01, cax=cax_01)
    cbar_01.ax.tick_params(labelsize=11)
    cbar_01.set_label('$|A|$', rotation=0, size=15, labelpad=-10, y=1.14)

    cax_04 = fig.add_axes([0.49, 0.11, 0.015, 0.37])
    cbar_04 = fig.colorbar(pc_04, cax=cax_04)
    cbar_04.ax.tick_params(labelsize=11)

    cax_03 = fig.add_axes([0.913, 0.11, 0.015, 0.37])
    cbar_03 = fig.colorbar(pc_03, cax=cax_03)
    cbar_03.ax.tick_params(labelsize=11)

    savedirectory = "C:/mnustes_science/simulation_data/FD/PDNLS_extended_PT/extras/dimensional/analysis/figures"
    plt.savefig('Fig01_04a.png', dpi=300)
    plt.close()