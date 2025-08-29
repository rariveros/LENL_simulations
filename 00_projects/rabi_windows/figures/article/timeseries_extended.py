from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    ######## ESTOS DATOS ESTAN EN EL DISCO EXTERNO #######
    disc = "D:/"
    #T = np.loadtxt(disc + r'mnustes_science\simulation_data\FD\rabi_windows\test\alpha=1.000\beta=1.000\mu=0.1000\nu=0.3200\sigma=3.000\gamma=0.2800\dist=20.0000\T.txt', delimiter=',')
    #X = np.loadtxt(disc + r'mnustes_science\simulation_data\FD\rabi_windows\test\alpha=1.000\beta=1.000\mu=0.1000\nu=0.3200\sigma=3.000\gamma=0.2800\dist=20.0000\X.txt', delimiter=',')
    #Zr = np.loadtxt(disc + r'mnustes_science\simulation_data\FD\rabi_windows\test\alpha=1.000\beta=1.000\mu=0.1000\nu=0.3200\sigma=3.000\gamma=0.2800\dist=20.0000\field_real.txt', delimiter=',')
    #Zi = np.loadtxt(disc + r'mnustes_science\simulation_data\FD\rabi_windows\test\alpha=1.000\beta=1.000\mu=0.1000\nu=0.3200\sigma=3.000\gamma=0.2800\dist=20.0000\field_img.txt', delimiter=',')

    #Z = Zr + 1j * Zi
    #Z_mod_ROs = np.abs(Z)

    T = np.loadtxt(disc + r'mnustes_science\simulation_data\FD\rabi_windows\test\alpha=1.000\beta=1.000\mu=0.1000\nu=0.3200\sigma=3.000\gamma=0.2800\dist=39.0000\T.txt', delimiter=',')
    X= np.loadtxt(disc + r'mnustes_science\simulation_data\FD\rabi_windows\test\alpha=1.000\beta=1.000\mu=0.1000\nu=0.3200\sigma=3.000\gamma=0.2800\dist=39.0000\X.txt', delimiter=',')
    Zr = np.loadtxt(disc + r'mnustes_science\simulation_data\FD\rabi_windows\test\alpha=1.000\beta=1.000\mu=0.1000\nu=0.3200\sigma=3.000\gamma=0.2800\dist=39.0000\field_real.txt', delimiter=',')
    Zi = np.loadtxt(disc + r'mnustes_science\simulation_data\FD\rabi_windows\test\alpha=1.000\beta=1.000\mu=0.1000\nu=0.3200\sigma=3.000\gamma=0.2800\dist=39.0000\field_img.txt', delimiter=',')

    Z = Zr + 1j * Zi

    Nx = len(X)
    Nt = len(T)
    Z_L = Z[:, 0:int(Nx / 2)]
    Z_R = Z[:, int(Nx / 2) + 1:]
    X_L = X[0:int(Nx / 2)]
    X_R = X[int(Nx / 2) + 1:]
    U_R = integrate.simpson(Z_R, X_R)
    U_L = integrate.simpson(Z_L, X_L)

    t0 = 2000
    ti, tf = 0, 400

    fig, (ax1, ax2) = plt.subplots(2,figsize=(3, 1.65))
    ax1.plot(T - t0, np.real(U_L), c="b", label="$\\textrm{Re }\Psi_L$", lw=2, zorder=5)
    ax1.plot(T - t0, np.imag(U_L), c="r", label="$\\textrm{Im }\Psi_L$", lw=2, zorder=5)
    ax1.plot(T - t0, np.abs(U_L), c="k", label="$|\Psi_L|$", lw=2, zorder=5)
    ax1.hlines(0, 0, 400, colors="k")
    ax1.set_ylabel('$\Psi_L$', size=15)
    ax1.set_xlim(0, 400)
    ax1.set_xticks([0, 200, 400])
    ax1.set_ylim(-0.6, 0.6)
    ax1.tick_params(axis="both", direction="in", labelsize=13, labelbottom=False)
    ax1.set_yticks([-0.5, 0, 0.5])
    #ax1.legend(fontsize=16)
    ax1.grid(alpha=0.3)

    ax2.plot(T - t0, np.real(U_R), c="b", label="$\\textrm{Re }\Psi_R$", lw=2, zorder=5)
    ax2.plot(T - t0, np.imag(U_R), c="r", label="$\\textrm{Im }\Psi_R$", lw=2, zorder=5)
    ax2.plot(T - t0, np.abs(U_R), c="k", label="$|\Psi_R|$", lw=2, zorder=5)
    ax2.hlines(0, 0, 400, colors="k")
    ax2.set_ylabel('$\Psi_R$', size=15)
    ax2.set_xlabel('$\\textrm{Time}$', size=15)
    ax2.set_xlim(0, 400)
    ax2.set_xticks([0, 200, 400])
    ax2.set_ylim(-0.6, 0.6)
    ax2.tick_params(axis="both", direction="in", labelsize=13)
    ax2.set_yticks([-0.5, 0, 0.5])
    ax2.grid(alpha=0.3)

    fig.subplots_adjust(left=0.23, right=0.95, bottom=0.3, top=0.95)
    plt.savefig("timeseries_extended.png", dpi=300)
    plt.close()