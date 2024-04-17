from back_process import *


if __name__ == "__main__":
    principal_sim = "C:/mnustes_science/simulation_data/FD/PDNLS_extended_PT/extras"
    principal_exp = "C:/mnustes_science/experimental_data"
    dir_02 = principal_sim + "/dimensional_dist_02/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.0180/sigma=6.000/gamma=0.1850/dist=31.000"#/dist_03/alpha=13.0480/beta=1.000/mu=0.0950/nu=0.0180/sigma=12.6700/gamma=0.1690/dist=40.000"#
    dir_01 = principal_sim + "/dimensional_dist_02/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.0180/sigma=6.000/gamma=0.1850/dist=21.0000"#"/dist_03/alpha=13.0480/beta=1.000/mu=0.0950/nu=0.0180/sigma=12.6700/gamma=0.1690/dist=31.000"#
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
    Z_modulo_01 = filtro_superficie(Z_modulo_01, 2, "XY") / np.sqrt(beta)

    ###
    xi_exp, xf_exp = -150, 150
    xi_sim, xf_sim = -50, 50
    ti_exp_01, tf_exp_01 = 0, 400
    ti_exp_02, tf_exp_02 = 0, 400
    ti_sim, tf_sim = 0, 400

    data = [[X_01, T_01 - (T_01[-1] - 400), Z_r_01], [X_01, T_01 - (T_01[-1] - 400), Z_i_01]]
    metadata = [[-150, 150], [-150, 150], [-100, 100], [-100, 100]]

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)

    pc_01 = ax1.pcolor(data[0][0], data[0][1], data[0][2], cmap=parula_map)
    ax1.set_ylim(0, 400)
    ax1.set_xlim(xi_sim, xf_sim)
    ax1.grid(alpha=0.2, color="k")
    ax1.tick_params(axis="x", direction="in", top=True, bottom=True)
    ax1.tick_params(axis="y", direction="in", left=True, right=True)
    ax1.fill_between([xi_sim, xi_sim + 0.225 * (xf_sim - xi_sim)], [0.705 * (tf_sim - ti_sim), tf_sim], tf_sim, color="white", zorder=12)
    ax1.plot([xi_sim, xi_sim + 0.23 * (xf_sim - xi_sim)], [0.7 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax1.plot([xi_sim, xi_sim + 0.23 * (xf_sim - xi_sim)], [tf_sim, tf_sim], linewidth=1, color="k", zorder=12)
    ax1.plot([xi_sim, xi_sim], [0.7 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax1.text(xi_sim + 0.005 * (xf_sim - xi_sim), ti_sim + 0.795 * (tf_sim - ti_sim), "$\\textbf{NUM}$", fontsize=13, rotation=45, zorder=13)
    ax1.text(xi_sim + 0.84 * (xf_sim - xi_sim), ti_sim + 0.87 * (tf_sim - ti_sim), "$\\textbf{(a)}$", fontsize=18, rotation=0, zorder=13, color="white")
    #ax1.text(xi_sim + 0.84 * (xf_sim - xi_sim), ti_sim + 0.87 * (tf_sim - ti_sim), "$\\textbf{(e)}$", fontsize=18, rotation=0, zorder=13, color="white")

    pc_02 = ax2.pcolor(data[1][0], data[1][1], data[1][2], cmap=parula_map)
    ax2.set_ylim(0, 400)
    ax2.set_yticklabels([])
    ax2.set_xlim(xi_sim, xf_sim)
    ax2.grid(alpha=0.2, color="k")
    ax2.tick_params(axis="x", direction="in", top=True, bottom=True)
    ax2.tick_params(axis="y", direction="in", left=True, right=True, zorder=10)
    ax2.fill_between([xi_sim, xi_sim + 0.225 * (xf_sim - xi_sim)], [0.705 * (tf_sim - ti_sim), tf_sim], tf_sim, color="white", zorder=12)
    ax2.plot([xi_sim, xi_sim + 0.23 * (xf_sim - xi_sim)], [0.7 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax2.plot([xi_sim, xi_sim + 0.23 * (xf_sim - xi_sim)], [tf_sim, tf_sim], linewidth=1, color="k", zorder=12)
    ax2.plot([xi_sim, xi_sim], [0.7 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax2.text(xi_sim + 0.005 * (xf_sim - xi_sim), ti_sim + 0.795 * (tf_sim - ti_sim), "$\\textbf{NUM}$", fontsize=13, rotation=45, zorder=13)
    ax2.text(xi_sim + 0.84 * (xf_sim - xi_sim), ti_sim + 0.87 * (tf_sim - ti_sim), "$\\textbf{(b)}$", fontsize=18, rotation=0, zorder=13, color="white")
    #ax2.text(xi_sim + 0.84 * (xf_sim - xi_sim), ti_sim + 0.87 * (tf_sim - ti_sim), "$\\textbf{(f)}$", fontsize=18, rotation=0, zorder=13, color="white")

    pc_03 = ax3.pcolor(data[0][0], data[0][1], data[0][2], cmap=parula_map)
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
    ax3.text(xi_sim + 0.84 * (xf_sim - xi_sim), ti_sim + 0.87 * (tf_sim - ti_sim), "$\\textbf{(a)}$", fontsize=18, rotation=0, zorder=13, color="white")
    #ax1.text(xi_sim + 0.84 * (xf_sim - xi_sim), ti_sim + 0.87 * (tf_sim - ti_sim), "$\\textbf{(e)}$", fontsize=18, rotation=0, zorder=13, color="white")

    pc_04 = ax4.pcolor(data[1][0], data[1][1], data[1][2], cmap=parula_map)
    ax4.set_ylim(0, 400)
    ax4.set_yticklabels([])
    ax4.set_xlim(xi_sim, xf_sim)
    ax4.grid(alpha=0.2, color="k")
    ax4.tick_params(axis="x", direction="in", top=True, bottom=True)
    ax4.tick_params(axis="y", direction="in", left=True, right=True, zorder=10)
    ax4.fill_between([xi_sim, xi_sim + 0.225 * (xf_sim - xi_sim)], [0.705 * (tf_sim - ti_sim), tf_sim], tf_sim, color="white", zorder=12)
    ax4.plot([xi_sim, xi_sim + 0.23 * (xf_sim - xi_sim)], [0.7 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax4.plot([xi_sim, xi_sim + 0.23 * (xf_sim - xi_sim)], [tf_sim, tf_sim], linewidth=1, color="k", zorder=12)
    ax4.plot([xi_sim, xi_sim], [0.7 * (tf_sim - ti_sim), tf_sim], linewidth=1, color="k", zorder=12)
    ax4.text(xi_sim + 0.005 * (xf_sim - xi_sim), ti_sim + 0.795 * (tf_sim - ti_sim), "$\\textbf{NUM}$", fontsize=13, rotation=45, zorder=13)
    ax4.text(xi_sim + 0.84 * (xf_sim - xi_sim), ti_sim + 0.87 * (tf_sim - ti_sim), "$\\textbf{(b)}$", fontsize=18, rotation=0, zorder=13, color="white")
    #ax2.text(xi_sim + 0.84 * (xf_sim - xi_sim), ti_sim + 0.87 * (tf_sim - ti_sim), "$\\textbf{(f)}$", fontsize=18, rotation=0, zorder=13, color="white")


    plt.subplots_adjust(wspace=0.2, hspace=0.07)
    #cax_01 = fig.add_axes([0.49, 0.51, 0.015, 0.37])
    #fig.colorbar(pc_01, cax=cax_01)
    cax_02 = fig.add_axes([0.913, 0.51, 0.015, 0.37])
    fig.colorbar(pc_02, cax=cax_02)
    savedirectory = "C:/mnustes_science/simulation_data/FD/PDNLS_extended_PT/extras/dimensional/analysis/figures"
    plt.savefig('Fig03a.png', dpi=300)
    plt.close()