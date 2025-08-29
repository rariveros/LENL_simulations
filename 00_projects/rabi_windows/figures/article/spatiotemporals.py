from back_process import *
from functions import *
from back_process import *
from time_integrators import *

if __name__ == "__main__":
    disc = "D:/"
    T = np.loadtxt(disc + r'mnustes_science\simulation_data\FD\rabi_windows\test\alpha=1.000\beta=1.000\mu=0.1000\nu=0.3200\sigma=3.000\gamma=0.2800\dist=20.0000\T.txt', delimiter=',')
    x_grid = np.loadtxt(disc + r'mnustes_science\simulation_data\FD\rabi_windows\test\alpha=1.000\beta=1.000\mu=0.1000\nu=0.3200\sigma=3.000\gamma=0.2800\dist=20.0000\X.txt', delimiter=',')
    Zr = np.loadtxt(disc + r'mnustes_science\simulation_data\FD\rabi_windows\test\alpha=1.000\beta=1.000\mu=0.1000\nu=0.3200\sigma=3.000\gamma=0.2800\dist=20.0000\field_real.txt', delimiter=',')
    Zi = np.loadtxt(disc + r'mnustes_science\simulation_data\FD\rabi_windows\test\alpha=1.000\beta=1.000\mu=0.1000\nu=0.3200\sigma=3.000\gamma=0.2800\dist=20.0000\field_img.txt', delimiter=',')

    Z = Zr + 1j * Zi
    Z_mod_ROs = np.abs(Z)

    T = np.loadtxt(disc + r'mnustes_science\simulation_data\FD\rabi_windows\test\alpha=1.000\beta=1.000\mu=0.1000\nu=0.3200\sigma=3.000\gamma=0.2800\dist=39.0000\T.txt', delimiter=',')
    x_grid = np.loadtxt(disc + r'mnustes_science\simulation_data\FD\rabi_windows\test\alpha=1.000\beta=1.000\mu=0.1000\nu=0.3200\sigma=3.000\gamma=0.2800\dist=39.0000\X.txt', delimiter=',')
    Zr = np.loadtxt(disc + r'mnustes_science\simulation_data\FD\rabi_windows\test\alpha=1.000\beta=1.000\mu=0.1000\nu=0.3200\sigma=3.000\gamma=0.2800\dist=39.0000\field_real.txt', delimiter=',')
    Zi = np.loadtxt(disc + r'mnustes_science\simulation_data\FD\rabi_windows\test\alpha=1.000\beta=1.000\mu=0.1000\nu=0.3200\sigma=3.000\gamma=0.2800\dist=39.0000\field_img.txt', delimiter=',')

    Z = Zr + 1j * Zi
    Z_mod_SP = np.abs(Z)

    t0 = 2000
    ti, tf = 0, 400
    xi, xf = -70, 70

    #################################
    fig, (ax1, ax2) = plt.subplots(figsize=(5, 5), nrows=2, ncols=1)

    pc_01 = ax1.pcolor(x_grid, T - t0, 100 * Z_mod_ROs, cmap="jet")
    ax1.set_ylim(ti, tf)
    ax1.set_xlim(xi, xf)
    ax1.set_yticks([0, 200, 400])
    ax1.set_ylabel("$\\textrm{Time}$", fontsize=20)
    ax1.grid(alpha=0.2, color="k")
    ax1.tick_params(axis="x", direction="in", labeltop=False, labelbottom=False, top=True, bottom=True, labelsize=15)
    ax1.tick_params(axis="y", direction="in", left=True, right=True, labelsize=15)

    pc_02 = ax2.pcolor(x_grid, T - t0, 100 * Z_mod_SP, cmap="jet")
    ax2.set_ylim(ti, tf)
    ax2.set_xlim(xi, xf)
    ax2.set_yticks([0, 200, 400])
    ax2.set_xlabel("$\\textrm{Space}$", fontsize=20)
    ax2.grid(alpha=0.2, color="k")
    ax2.tick_params(axis="x", direction="in", labeltop=False, labelbottom=True, top=True, bottom=True, labelsize=15)
    ax2.tick_params(axis="y", direction="in", left=True, right=True, labelleft=True, labelsize=15)

    plt.subplots_adjust(wspace=0.2, hspace=0.15, left=0.23, right=0.8, bottom=0.2, top=0.95)

    cax_01 = fig.add_axes([0.82, 0.6, 0.025, 0.35])
    fig.colorbar(pc_01, cax=cax_01)
    cb1 = fig.colorbar(pc_01, cax=cax_01)
    cb1.ax.tick_params(labelsize=14)

    cax_02 = fig.add_axes([0.82, 0.2, 0.025, 0.35])
    fig.colorbar(pc_02, cax=cax_02)
    cb2 = fig.colorbar(pc_02, cax=cax_02)
    cb2.ax.tick_params(labelsize=14)

    plt.savefig("ST.png", dpi=300)
    plt.close()