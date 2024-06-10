import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    disco = 'D:/'
    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FDa'
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    directories = [name for name in os.listdir(working_directory) if os.path.isdir(os.path.join(working_directory, name))]

    datafile = "C:/mnustes_science/simulation_data/FD/PDNLS_extended_PT/extras/dimensional/analysis"
    fig, ((ax1, ax2, ax3)) = plt.subplots(nrows=3, ncols=1)
    if not os.path.exists(datafile):
        os.makedirs(datafile)
    ZR = []
    ZI = []
    ZR_ansatz = []
    ZI_ansatz = []
    gammas = []
    n = 0
    for directory in directories:
        print("#############   " + directory + "   #############")
        data = np.loadtxt(working_directory + "/" + directory + '/data.txt', delimiter=',')
        params = np.loadtxt(working_directory + "/" + directory + '/parameters.txt', delimiter=',')
        Z_r = np.loadtxt(working_directory + "/" + directory + '/field_real.txt', delimiter=',')
        Z_i = np.loadtxt(working_directory + "/" + directory + '/field_img.txt', delimiter=',')
        X = np.loadtxt(working_directory + "/" + directory + '/X.txt', delimiter=',')
        #[alpha, beta, gamma_0, mu_0, nu, sigma]
        #[center, delta, ancho, amplitud]
        center = data[0]
        centerR = data[1]
        centerI = data[2]

        alpha = params[0]
        gamma_0 = params[2]
        mu = params[3]
        nu = params[4]
        sigma = params[5]

        gammaR_ij = gamma_0 * (np.exp(-centerR ** 2 / (2 * sigma ** 2)))
        gammaI_ij = gamma_0 * (np.exp(-centerI ** 2 / (2 * sigma ** 2)))
        phaseR_ij = 0.5 * np.arccos(mu / gammaR_ij)
        phaseI_ij = 0.5 * np.arccos(mu / gammaI_ij)
        deltaR_ij = (- nu + np.sqrt(gammaR_ij ** 2 - mu ** 2))
        deltaI_ij = 0.77 * (- nu + np.sqrt(gammaI_ij ** 2 - mu ** 2))
        AR_i = np.sqrt(2 * deltaR_ij) * np.cos(phaseR_ij)
        AI_i = - np.sqrt(2 * deltaI_ij) * np.sin(phaseI_ij)
        phi_01_i = (AR_i / np.cosh(np.sqrt(deltaR_ij / alpha) * (X - centerR)))
        phi_02_i = (AI_i / np.cosh(np.sqrt(deltaI_ij / alpha) * (X - centerI)))
        ax1.plot(X, Z_r[-1, :] + n, c="b")
        ax1.plot(X, phi_01_i + n, c="b", linestyle="--")
        ax2.plot(X, Z_i[-1, :] + n, c="r")
        ax2.plot(X, phi_02_i + n, c="r", linestyle="--")
        ax3.plot(X, np.abs(Z_r[-1, :] + 1j * Z_i[-1, :]) + n, c="k")
        ax3.plot(X, np.abs(phi_01_i + 1j * phi_02_i) + n, c="k", linestyle="--")
        ZR.append(Z_r[-1, :])
        ZI.append(Z_i[-1, :])
        ZR_ansatz.append(phi_01_i)
        ZI_ansatz.append(phi_02_i)
        gammas.append(gamma_0)
        n = n + 0.4
    plt.show()
    plt.close()

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)

    pc_01 = ax1.pcolor(X, gammas, ZR, cmap='nipy_spectral')
    ax1.set_ylabel("$\gamma_0$", fontsize=15)
    ax1.tick_params(axis="x", direction="in", labeltop=True, labelbottom=False, top=True, bottom=True)
    ax1.tick_params(axis="y", direction="in", left=True, right=True)
    ax1.grid(alpha=0.2, color="k")

    pc_02 = ax2.pcolor(X, gammas, np.abs(ZI), cmap='nipy_spectral')
    ax2.set_yticklabels([])
    ax2.tick_params(axis="x", direction="in", labeltop=True, labelbottom=False, top=True, bottom=True)
    ax2.tick_params(axis="y", direction="in", left=True, right=True)
    ax2.grid(alpha=0.2, color="k")

    pc_03 = ax3.pcolor(X, gammas, ZR_ansatz, cmap='nipy_spectral')
    ax3.set_ylabel("$\gamma_0$", fontsize=15)
    ax3.tick_params(axis="x", direction="in", labeltop=False, labelbottom=True, top=True, bottom=True)
    ax3.tick_params(axis="y", direction="in", left=True, right=True)
    ax3.grid(alpha=0.2, color="k")

    pc_04 = ax4.pcolor(X, gammas, np.abs(ZI_ansatz), cmap='nipy_spectral')
    ax4.set_yticklabels([])
    ax4.tick_params(axis="x", direction="in", labeltop=False, labelbottom=True, top=True, bottom=True)
    ax4.tick_params(axis="y", direction="in", left=True, right=True)
    ax4.grid(alpha=0.2, color="k")

    plt.subplots_adjust(wspace=0.25, hspace=0.07)
    cax_01 = fig.add_axes([0.49, 0.51, 0.015, 0.37])
    cbar_01 = fig.colorbar(pc_01, cax=cax_01)
    cbar_01.ax.tick_params(labelsize=11)

    cax_02 = fig.add_axes([0.913, 0.51, 0.015, 0.37])
    cbar_02 = fig.colorbar(pc_02, cax=cax_02)
    cbar_02.ax.tick_params(labelsize=11)

    cax_03 = fig.add_axes([0.49, 0.11, 0.015, 0.37])
    cbar_03 = fig.colorbar(pc_03, cax=cax_03)
    cbar_03.ax.tick_params(labelsize=11)

    cax_04 = fig.add_axes([0.913, 0.11, 0.015, 0.37])
    cbar_04 = fig.colorbar(pc_04, cax=cax_04)
    cbar_04.ax.tick_params(labelsize=11)

    plt.show()
    plt.close()
