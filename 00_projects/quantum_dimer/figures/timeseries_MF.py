from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    ######## ESTOS DATOS ESTAN EN EL DISCO EXTERNO #######
    #directory = r"C:\mnustes_science\simulation_data\FD\coherent\meanfield\Delta=0.0000\gamma=0.1000\Omega=0.0200\k=0.3000" #DAMPED
    #directory = r"C:\mnustes_science\simulation_data\FD\coherent\meanfield\Delta=0.0000\gamma=0.1000\Omega=0.1000\k=0.3000"  # ROs
    #directory = r"C:\mnustes_science\simulation_data\FD\coherent\meanfield\Delta=0.0000\gamma=0.1000\Omega=0.1500\k=0.3000"  # Strange_ROs
    directory = r"C:\mnustes_science\simulation_data\FD\coherent\meanfield\Delta=0.0000\gamma=0.1000\Omega=0.2400\k=0.3000"  # SP

    T = np.loadtxt(directory + r'\T.txt', delimiter=',')[:-1]
    U = np.loadtxt(directory + r'\U.txt', delimiter=',', dtype=complex)
    V = np.loadtxt(directory + r'\V.txt', delimiter=',', dtype=complex)

    amp_factor = 10

    t0 = 200
    ti, tf = 0, 200

    deep_blue = [1 / 255, 13 / 255, 124 / 255]
    deep_red= [173 / 255, 0, 0]

    fig, (ax1, ax2) = plt.subplots(2,figsize=(3, 1.65))
    ax1.plot(T - t0, np.real(U) * amp_factor, c=deep_blue, label=r"$\textrm{Re }\alpha_1$", lw=2, zorder=5)
    ax1.plot(T - t0, np.imag(U) * amp_factor, c=deep_red, label=r"$\textrm{Im }\alpha_1$", lw=2, zorder=5)
    #ax1.plot(T - t0, np.abs(U), c="k", label=r"$|\Psi_L|$", lw=2, zorder=5)
    ax1.hlines(0, 0, 400, colors="k")
    ax1.set_ylabel(r'$\alpha_1$', size=15)
    ax1.set_xlim(0, tf)
    ax1.set_xticks([0, tf])
    ax1.set_ylim(-6.5, 6.5)
    ax1.tick_params(axis="both", direction="in", labelsize=13, labelbottom=False)
    ax1.set_yticks([-4, 0, 4])
    #ax1.legend(fontsize=16)
    ax1.grid(alpha=0.3)

    ax2.plot(T - t0, np.real(V) * amp_factor, c=deep_blue, label=r"$\textrm{Re }\alpha_2$", lw=2, zorder=5)
    ax2.plot(T - t0, np.imag(V) * amp_factor, c=deep_red, label=r"$\textrm{Im }\alpha_2$", lw=2, zorder=5)
    #ax2.plot(T - t0, np.abs(V), c="k", label="$|\Psi_R|$", lw=2, zorder=5)
    ax2.hlines(0, 0, 400, colors="k")
    ax2.set_ylabel(r'$\alpha_2$', size=15)
    ax2.set_xlabel(r'$\textrm{Time}$', size=15)
    ax2.set_xlim(0, tf)
    ax2.set_ylim(-6.5, 6.5)
    ax2.set_xticks([0, tf])
    ax2.set_yticks([-4, 0, 4])
    ax2.tick_params(axis="both", direction="in", labelsize=13)
    ax2.grid(alpha=0.3)

    fig.subplots_adjust(left=0.23, right=0.95, bottom=0.3, top=0.95)
    plt.savefig("timeseries.png", dpi=300)
    plt.close()