from functions import *
from back_process import *
import matplotlib.colors as mcolors

if __name__ == '__main__':
    Ns = np.arange(0, 20.1, 0.1)

    k = 1
    k1 = k
    gamma1 = np.arange(-0.99, 0.99, 0.01) * k
    N1, GAMMA1 = np.meshgrid(Ns, gamma1)
    F1a = 2 * np.sqrt(k ** 2 - GAMMA1 ** 2 + 0j) * (
                (2 * 10 * N1) / ((N1 + 2) ** 2 * np.sqrt(k ** 2 - GAMMA1 ** 2 + 0j)) - 1) ** 0.5
    F1a_r = np.real(F1a);
    F1a_i = np.imag(F1a);
    F1a_plot = F1a_r - F1a_i
    F1b = 2 * np.sqrt(k ** 2 - GAMMA1 ** 2 + 0j) * (
                -(2 * 10 * N1) / ((N1 + 2) ** 2 * np.sqrt(k ** 2 - GAMMA1 ** 2 + 0j)) - 1) ** 0.5
    F1b_r = np.real(F1b);
    F1b_i = np.imag(F1b);
    F1b_plot = F1b_r - F1b_i

    k = 2
    k2 = k
    gamma2 = np.arange(-0.99, 0.99, 0.01) * k
    N2, GAMMA2 = np.meshgrid(Ns, gamma2)
    F2a = 2 * np.sqrt(k ** 2 - GAMMA2 ** 2 + 0j) * (
                (2 * 10 * N2) / ((N2 + 2) ** 2 * np.sqrt(k ** 2 - GAMMA2 ** 2 + 0j)) - 1) ** 0.5
    F2a_r = np.real(F2a);
    F2a_i = np.imag(F2a);
    F2a_plot = F2a_r - F2a_i
    F2b = 2 * np.sqrt(k ** 2 - GAMMA2 ** 2 + 0j) * (
                -(2 * 10 * N2) / ((N2 + 2) ** 2 * np.sqrt(k ** 2 - GAMMA2 ** 2 + 0j)) - 1) ** 0.5
    F2b_r = np.real(F2b);
    F2b_i = np.imag(F2b);
    F2b_plot = F2b_r - F2b_i

    k = 2.5
    k3 = k
    gamma3 = np.arange(-0.99, 0.99, 0.01) * k
    N3, GAMMA3 = np.meshgrid(Ns, gamma3)
    F3a = 2 * np.sqrt(k ** 2 - GAMMA3 ** 2 + 0j) * (
                (2 * 10 * N3) / ((N3 + 2) ** 2 * np.sqrt(k ** 2 - GAMMA3 ** 2 + 0j)) - 1) ** 0.5
    F3a_r = np.real(F3a);
    F3a_i = np.imag(F3a);
    F3a_plot = F3a_r - F3a_i
    F3b = 2 * np.sqrt(k ** 2 - GAMMA3 ** 2 + 0j) * (
                -(2 * 10 * N3) / ((N3 + 2) ** 2 * np.sqrt(k ** 2 - GAMMA3 ** 2 + 0j)) - 1) ** 0.5
    F3b_r = np.real(F3b);
    F3b_i = np.imag(F3b);
    F3b_plot = F3b_r - F3b_i

    fig, axes = plt.subplots(3, 2, figsize=(6, 6))

    labelsize = 15
    xy_labelsize = 18

    ################## k = 1.0 ##################

    F1a_plot[0, 0] = -np.abs(np.amax(F1a_plot));
    F1b_plot[0, 0] = np.abs(np.amin(F1b_plot))
    pc_01a = axes[0, 0].pcolor(GAMMA1/k1, N1, -F1a_plot, cmap='bwr')
    axes[0, 0].tick_params(axis="y", direction="in", labelsize=labelsize, left=True, right=True, labelleft=True,
                           labelright=False)
    axes[0, 0].tick_params(axis="x", direction="in", labelsize=labelsize, top=True, bottom=True, labeltop=False,
                           labelbottom=False)
    pc_01b = axes[0, 1].pcolor(GAMMA1/k1, N1, -F1b_plot, cmap='bwr')
    axes[0, 1].tick_params(axis="y", direction="in", labelsize=labelsize, left=True, right=True, labelleft=False,
                           labelright=False)
    axes[0, 1].tick_params(axis="x", direction="in", labelsize=labelsize, top=True, bottom=True, labeltop=False,
                           labelbottom=False)
    axes[0, 0].set_ylabel("$N_{*}$", fontsize=xy_labelsize)
    axes[0, 0].set_xticks([-1, 0, 1])
    axes[0, 0].set_yticks([0, 10, 20])
    axes[0, 1].set_xticks([-1, 0, 1])
    axes[0, 1].set_yticks([0, 10, 20])

    cax_01a = fig.add_axes([0.45, 0.65, 0.015, 0.2]);
    cb1a = fig.colorbar(pc_01a, cax=cax_01a);
    cb1a.ax.tick_params(labelsize=labelsize)
    cax_01b = fig.add_axes([0.91, 0.65, 0.015, 0.2]);
    cb1b = fig.colorbar(pc_01b, cax=cax_01b);
    cb1b.ax.tick_params(labelsize=labelsize)
    cb1a.set_ticks([-2, 0, 2])  # choose tick positions
    cb1a.set_ticklabels(["$2$", "$0$", "$2i$"])  # custom names
    cb1b.set_ticks([-3, 0, 3])  # choose tick positions
    cb1b.set_ticklabels(["$3$", "$0$", "$3i$"])  # custom names
    cb1a.set_label('$\lambda_{+}^{(1)}$', rotation=0, size=xy_labelsize, labelpad=-18, y=1.32)
    cb1b.set_label('$\lambda_{-}^{(1)}$', rotation=0, size=xy_labelsize, labelpad=-18, y=1.32)

    ################## k = 2.0 ##################

    F2a_plot[0, 0] = np.abs(np.amin(F2a_plot));
    F2b_plot[0, 0] = np.abs(np.amin(F2b_plot))
    pc_02a = axes[1, 0].pcolor(GAMMA2/k2, N2, -F2a_plot, cmap='bwr')
    axes[1, 0].tick_params(axis="y", direction="in", labelsize=labelsize, left=True, right=True, labelleft=True,
                           labelright=False)
    axes[1, 0].tick_params(axis="x", direction="in", labelsize=labelsize, top=True, bottom=True, labeltop=False,
                           labelbottom=False)
    pc_02b = axes[1, 1].pcolor(GAMMA2/k2, N2, -F2b_plot, cmap='bwr')
    axes[1, 1].tick_params(axis="y", direction="in", labelsize=labelsize, left=True, right=True, labelleft=False,
                           labelright=False)
    axes[1, 1].tick_params(axis="x", direction="in", labelsize=labelsize, top=True, bottom=True, labeltop=False,
                           labelbottom=False)
    axes[1, 0].set_ylabel("$N_{*}$", fontsize=xy_labelsize)
    axes[1, 0].set_xticks([-1, 0, 1])
    axes[1, 0].set_yticks([0, 10, 20])
    axes[1, 1].set_xticks([-1, 0, 1])
    axes[1, 1].set_yticks([0, 10, 20])
    cax_02a = fig.add_axes([0.45, 0.40, 0.015, 0.2]);
    cb2a = fig.colorbar(pc_02a, cax=cax_02a);
    cb2a.ax.tick_params(labelsize=labelsize)
    cax_02b = fig.add_axes([0.91, 0.4, 0.015, 0.2]);
    cb2b = fig.colorbar(pc_02b, cax=cax_02b);
    cb2b.ax.tick_params(labelsize=labelsize)
    cb2a.set_ticks([-3, 0, 3])  # choose tick positions
    cb2a.set_ticklabels(["$3$", "$0$", "$3i$"])  # custom names
    cb2b.set_ticks([-5, 0, 5])  # choose tick positions
    cb2b.set_ticklabels(["$5$", "$0$", "$5i$"])  # custom names

    ################## k = 2.5 ##################

    F3a_plot[0, 0] = np.abs(np.amin(F3a_plot));
    F3b_plot[0, 0] = np.abs(np.amin(F3b_plot))
    pc_03a = axes[2, 0].pcolor(GAMMA3/k3, N3, -F3a_plot, cmap='bwr')
    axes[2, 0].tick_params(axis="y", direction="in", labelsize=labelsize, left=True, right=True, labelleft=True,
                           labelright=False)
    axes[2, 0].tick_params(axis="x", direction="in", labelsize=labelsize, top=True, bottom=True, labeltop=False,
                           labelbottom=True)
    pc_03b = axes[2, 1].pcolor(GAMMA3/k3, N3, -F3b_plot, cmap='bwr')
    axes[2, 1].tick_params(axis="y", direction="in", labelsize=labelsize, left=True, right=True, labelleft=False,
                           labelright=False)
    axes[2, 1].tick_params(axis="x", direction="in", labelsize=labelsize, top=True, bottom=True, labeltop=False,
                           labelbottom=True)
    axes[2, 0].set_xlabel("$\gamma/\kappa$", fontsize=xy_labelsize)
    axes[2, 0].set_ylabel("$N_{*}$", fontsize=xy_labelsize)
    axes[2, 1].set_xlabel("$\gamma/\kappa$", fontsize=xy_labelsize)
    axes[2, 0].set_xticks([-1, 0, 1])
    axes[2, 0].set_yticks([0, 10, 20])
    axes[2, 1].set_xticks([-1, 0, 1])
    axes[2, 1].set_yticks([0, 10, 20])
    cax_03a = fig.add_axes([0.45, 0.15, 0.015, 0.2]);
    cb3a = fig.colorbar(pc_03a, cax=cax_03a);
    cb3a.ax.tick_params(labelsize=labelsize)
    cax_03b = fig.add_axes([0.91, 0.145, 0.015, 0.2]);
    cb3b = fig.colorbar(pc_03b, cax=cax_03b);
    cb3b.ax.tick_params(labelsize=labelsize)
    cb3a.set_ticks([-5, 0, 5])  # choose tick positions
    cb3a.set_ticklabels(["$5$", "$0$", "$5i$"])  # custom names
    cb3b.set_ticks([-6, 0, 6])  # choose tick positions
    cb3b.set_ticklabels(["$6$", "$0$", "$6i$"])  # custom names

    fig.subplots_adjust(wspace=0.35, hspace=0.25, left=0.1, right=0.9, bottom=0.15, top=0.85)
    plt.savefig("eigen_shit.png", dpi=300)
    plt.close()