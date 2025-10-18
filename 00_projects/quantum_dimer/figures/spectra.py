from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    working_directory = r"C:\mnustes_science\simulation_data\FD\coherent\spectrums_01\Delta=0.0000\gamma=0.1000\Omega=0.2000"
    directories = [name for name in os.listdir(working_directory) if os.path.isdir(os.path.join(working_directory, name))]
    dir_0 = working_directory + "/" + directories[0]
    dir_1 = working_directory + "/" + directories[1]
    dir_2 = working_directory + "/" + directories[2]
    dir_3 = working_directory + "/" + directories[3]

    S0 = np.loadtxt(dir_0 + "/Spectra_V.txt", delimiter=',')
    print("### Spectrum #0 Ready! ###")
    S1 = np.loadtxt(dir_1 + "/Spectra_V.txt", delimiter=',')
    print("### Spectrum #1 Ready! ###")
    S2 = np.loadtxt(dir_2 + "/Spectra_V.txt", delimiter=',')
    print("### Spectrum #2 Ready! ###")
    S3 = np.loadtxt(dir_3 + "/Spectra_V.txt", delimiter=',')
    print("### Spectrum #3 Ready! ###")

    S0 = filtro_superficie(S0, 2, "Y")
    S0 = filtro_superficie(S0, 2, "X")
    S1 = filtro_superficie(S1, 2, "Y")
    S1 = filtro_superficie(S1, 2, "X")
    S2 = filtro_superficie(S2, 2, "Y")
    S2 = filtro_superficie(S2, 2, "X")
    S3 = filtro_superficie(S3, 2, "Y")
    S3 = filtro_superficie(S3, 2, "X")


    F0 = np.loadtxt(dir_0 + "/Freqs.txt", delimiter=',')
    F1 = np.loadtxt(dir_1 + "/Freqs.txt", delimiter=',')
    F2 = np.loadtxt(dir_2 + "/Freqs.txt", delimiter=',')
    F3 = np.loadtxt(dir_3 + "/Freqs.txt", delimiter=',')

    K0 = np.loadtxt(dir_0 + "/Ks_mean.txt", delimiter=',')
    K1 = np.loadtxt(dir_1 + "/Ks_mean.txt", delimiter=',')
    K2 = np.loadtxt(dir_2 + "/Ks_mean.txt", delimiter=',')
    K3 = np.loadtxt(dir_3 + "/Ks_mean.txt", delimiter=',')

    fig, (ax0, ax1, ax2, ax3) = plt.subplots(4, 1, figsize=(4, 6.5))
    ax0.pcolormesh(K0, F0, S0.T,shading='auto', cmap='inferno', norm=LogNorm(vmin=1e4, vmax=1e12))
    ax0.set_ylabel(r"$\omega$", fontsize=15)
    ax0.set_yticks([-0.3, 0, 0.3])
    ax0.tick_params(axis="both", direction="in", labelsize=11, top=True, bottom=True, labeltop=False, labelbottom=False)
    ax0.set_ylim(-0.4, 0.4)

    ax1.pcolormesh(K1, F1, S1.T,shading='auto', cmap='inferno', norm=LogNorm(vmin=1e4, vmax=1e12))
    ax1.set_ylabel(r"$\omega$", fontsize=15)
    ax1.set_yticks([-0.3, 0, 0.3])
    ax1.tick_params(axis="both", direction="in", labelsize=11, top=True, bottom=True, labeltop=False, labelbottom=False)
    ax1.set_ylim(-0.4, 0.4)

    ax2.pcolormesh(K2, F2, S2.T,shading='auto', cmap='inferno', norm=LogNorm(vmin=1e4, vmax=1e12))
    ax2.set_ylabel(r"$\omega$", fontsize=15)
    ax2.set_yticks([-0.3, 0, 0.3])
    ax2.tick_params(axis="both", direction="in", labelsize=11, top=True, bottom=True, labeltop=False, labelbottom=False)
    ax2.set_ylim(-0.4, 0.4)

    pcm = ax3.pcolormesh(K3, F3, S3.T,shading='auto', cmap='inferno', norm=LogNorm(vmin=1e4, vmax=1e12))
    ax3.set_ylabel(r"$\omega$", fontsize=15)
    ax3.set_xlabel(r"$\kappa$", fontsize=15)
    ax3.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax3.set_yticks([-0.3, 0, 0.3])
    ax3.tick_params(axis="both", direction="in", labelsize=11, top=True, bottom=True, labeltop=False, labelbottom=True)
    ax3.set_ylim(-0.4, 0.4)

    cax = fig.add_axes([0.82, 0.11, 0.025, 0.77])
    cb = fig.colorbar(pcm, cax=cax)
    cb.set_label(r"$\textrm{log}_{10}\ S_1(\omega)$", fontsize=14)
    cb.ax.tick_params(labelsize=10)
    fig.subplots_adjust(wspace=0.28, hspace=0.1, left=0.15, right=0.8) #left=0.1, right=0.9, bottom=0.25, top=0.8)

    plt.savefig("spectra.png", dpi=300)
    plt.close()