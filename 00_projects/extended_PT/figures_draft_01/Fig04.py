from back_process import *


if __name__ == "__main__":
    principal_sim = "C:/mnustes_science/simulation_data/FD/PDNLS_extended_PT/extras"
    principal_exp = "C:/mnustes_science/experimental_data"
    dir_01 = principal_sim + "/dimensional/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.0180/sigma=6.000/gamma=0.2800/dist=20.000"
    dir_02 = principal_sim + "/dimensional/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.0180/sigma=6.000/gamma=0.2200/dist=20.000"
    dir_03 = principal_exp + "/PT_05/f=14.20_a=12.00"
    dir_04 = principal_exp + "/PT_05/f=14.20_a=14.00"

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

    x03 = injection_03[0] + injection_03[1] - 16 * 2.5
    X_03 = X_03 - x03
    T_03 = T_03 * (w_1 / (2 * np.pi))

    Z_04 = np.loadtxt(dir_04 + '/Z_strobo.txt', delimiter=',')
    X_04 = np.loadtxt(dir_04 + '/X_mm.txt', delimiter=',')
    T_04 = np.loadtxt(dir_04 + '/T_strobo.txt', delimiter=',')
    injection_04 = np.loadtxt(dir_04 + '/IL_mm.txt', delimiter=',')
    Z_04 = filtro_superficie(Z_04, 5, "YX")

    x04 = injection_04[0] + injection_04[1] - 16 * 3
    X_04 = X_04 - x04
    T_04 = T_04 * (w_1 / (2 * np.pi))

    xi_exp, xf_exp = -150, 150
    xi_sim, xf_sim = -90, 90
    ti_exp_01, tf_exp_01 = 0, 200
    ti_exp_02, tf_exp_02 = 0, 350
    ti_sim, tf_sim = 0, 400

    data = [[X_04, T_04 - (T_04[-1] - 400), Z_04], [X_03, T_03 - (T_03[-1] - 400), Z_03], [X_01, T_01 - (T_01[-1] - 400), np.flip(Z_modulo_01, axis=1)], [X_02, T_02 - (T_02[-1] - 400), Z_modulo_02]]
    metadata = [[-150, 150], [-150, 150], [-100, 100], [-100, 100]]

    fig, ((ax2, ax1), (ax4, ax3)) = plt.subplots(nrows=2, ncols=2)

    ax1.pcolor(data[0][0], data[0][1], data[0][2], cmap=parula_map)
    ax1.set_ylim(ti_exp_01, tf_exp_01)
    ax1.set_xlim(-150, 150)
    ax1.set_yticklabels([])
    ax1.tick_params(axis="x", direction="in", labeltop=True, labelbottom=False, top=True, bottom=True)
    ax1.tick_params(axis="y", direction="in", left=True, right=True)
    ax1.grid(alpha=0.2, color="k")
    ax1.fill_between([xi_exp, xi_exp + 0.225 * (xf_exp - xi_exp)], [0.705 * (tf_exp_01 - ti_exp_01), tf_exp_01], tf_exp_01, color="white", zorder=12)
    ax1.plot([xi_exp, xi_exp + 0.23 * (xf_exp - xi_exp)], [0.7 * (tf_exp_01 - ti_exp_01), tf_exp_01], linewidth=1, color="k", zorder=12)
    ax1.plot([xi_exp, xi_exp + 0.23 * (xf_exp - xi_exp)], [tf_exp_01, tf_exp_01], linewidth=1, color="k", zorder=12)
    ax1.plot([xi_exp, xi_exp], [0.7 * (tf_exp_01 - ti_exp_01), tf_exp_01], linewidth=1, color="k", zorder=12)
    ax1.text(xi_exp + 0.005 * (xf_exp - xi_exp), ti_exp_01 + 0.795 * (tf_exp_01 - ti_exp_01), "$\\textbf{EXP}$", fontsize=13, rotation=45, zorder=13)

    ax2.pcolor(data[1][0], data[1][1], data[1][2], cmap=parula_map)
    ax2.set_ylim(ti_exp_02, tf_exp_02)

    ax2.set_xlim(-150, 150)
    ax2.tick_params(axis="x", direction="in", labeltop=True, labelbottom=False, top=True, bottom=True)
    ax2.tick_params(axis="y", direction="in", left=True, right=True)
    ax2.grid(alpha=0.2, color="k")
    ax2.fill_between([xi_exp, xi_exp + 0.225 * (xf_exp - xi_exp)], [0.705 * (tf_exp_02 - ti_exp_02), tf_exp_02], tf_exp_02, color="white", zorder=12)
    ax2.plot([xi_exp, xi_exp + 0.23 * (xf_exp - xi_exp)], [0.7 * (tf_exp_02 - ti_exp_02), tf_exp_02], linewidth=1, color="k", zorder=12)
    ax2.plot([xi_exp, xi_exp + 0.23 * (xf_exp - xi_exp)], [tf_exp_02, tf_exp_02], linewidth=1, color="k", zorder=12)
    ax2.plot([xi_exp, xi_exp], [0.7 * (tf_exp_02 - ti_exp_02), tf_exp_02], linewidth=1, color="k", zorder=12)
    ax2.text(xi_exp + 0.005 * (xf_exp - xi_exp), ti_exp_02 + 0.795 * (tf_exp_02 - ti_exp_02), "$\\textbf{EXP}$", fontsize=13, rotation=45, zorder=13)

    ax3.pcolor(data[2][0], data[2][1], data[2][2], cmap=parula_map)
    ax3.set_ylim(0, 400)
    ax3.set_xlim(-90, 90)
    ax3.grid(alpha=0.2, color="k")
    ax3.tick_params(axis="x", direction="in", top=True, bottom=True)
    ax3.tick_params(axis="y", direction="in", left=True, right=True)
    ax3.set_yticklabels([])
    ax3.fill_between([xi_sim, xi_sim + 0.225 * (xf_sim - xi_sim)], [0.705 * (tf_sim - ti_sim), tf_sim], tf_sim, color="white", zorder=12)
    ax3.plot([xi_sim, xi_sim + 0.23 * (xf_sim - xi_sim)], [0.7 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax3.plot([xi_sim, xi_sim + 0.23 * (xf_sim - xi_sim)], [tf_sim, tf_sim], linewidth=1, color="k", zorder=12)
    ax3.plot([xi_sim, xi_sim], [0.7 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax3.text(xi_sim + 0.005 * (xf_sim - xi_sim), ti_sim + 0.795 * (tf_sim - ti_sim), "$\\textbf{NUM}$", fontsize=13, rotation=45, zorder=13)

    pc = ax4.pcolor(data[3][0], data[3][1], data[3][2], cmap=parula_map)
    ax4.set_ylim(0, 400)
    ax4.set_xlim(-90, 90)
    ax4.grid(alpha=0.1, color="k")
    ax4.tick_params(axis="x", direction="in", top=True, bottom=True)
    ax4.tick_params(axis="y", direction="in", left=True, right=True, zorder=10)
    ax4.fill_between([xi_sim, xi_sim + 0.225 * (xf_sim - xi_sim)], [0.705 * (tf_sim - ti_sim), tf_sim], tf_sim, color="white", zorder=12)
    ax4.plot([xi_sim, xi_sim + 0.23 * (xf_sim - xi_sim)], [0.7 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax4.plot([xi_sim, xi_sim + 0.23 * (xf_sim - xi_sim)], [tf_sim, tf_sim], linewidth=1, color="k", zorder=12)
    ax4.plot([xi_sim, xi_sim], [0.7 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax4.text(xi_sim + 0.005 * (xf_sim - xi_sim), ti_sim + 0.795 * (tf_sim - ti_sim), "$\\textbf{NUM}$", fontsize=13, rotation=45, zorder=13)

    plt.subplots_adjust(wspace=0.07, hspace=0.07)
    cax = fig.add_axes([0.92, 0.11, 0.015, 0.77])
    fig.colorbar(pc, cax=cax)
0.