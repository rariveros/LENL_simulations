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
    centers = []
    nus = []
    gammas = ["0.180"]
    sigmas = ["15.000", "17.500", "20.000", "22.500", "25.000", "27.500", "30.000"]
    colors = (np.array([10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30]) - 10) / 20
    fig, ((ax1, ax2),(ax3, ax4)) = plt.subplots(nrows=2, ncols=2)
    if not os.path.exists(datafile):
        os.makedirs(datafile)
    for gamma in gammas:
        for i in range(len(sigmas)):
            sigma = sigmas[i]
            deltas = []
            amps = []
            anchos = []
            for directory in directories:
                if os.path.exists(working_directory + "/" + directory + '/sigma=' + sigma + '/gamma=' + gamma):
                    print("#############   " + directory + "   #############")
                    data = np.loadtxt(working_directory + "/" + directory + '/sigma=' + sigma + '/gamma=' + gamma + '/data.txt', delimiter=',')
                    params = np.loadtxt(working_directory + "/" + directory + '/sigma=' + sigma + '/gamma=' + gamma + '/parameters.txt', delimiter=',')
                    #[alpha, beta, gamma_0, mu_0, nu, sigma]
                    #[center, delta, ancho, amplitud]
                    delta = data[1]
                    amp = data[3]
                    ancho = data[2]
                    deltas.append(np.abs(delta))
                    amps.append(np.abs(amp))
                    anchos.append(np.abs(ancho))
                    color = (colors[i], 0, 1 - colors[i])
                    depth = 3
                    edge = "black"
                    if len(data) == 5:
                        color = (1, 1, 1)
                        depth = 1
                        edge = "white"
                    if params[4] == -0.15:
                        ax1.scatter(np.abs(amp), np.abs(delta), label="$\sigma_i=" + sigma + "$", zorder=depth, c=color, edgecolors=edge, alpha=1)
                        ax2.scatter(np.abs(ancho), np.abs(delta), label="$\sigma_i=" + sigma + "$", zorder=depth, c=color, edgecolors=edge, alpha=1)
                    else:
                        ax1.scatter(np.abs(amp), np.abs(delta), zorder=depth, c=color, edgecolors=edge, alpha=1)
                        ax2.scatter(np.abs(ancho), np.abs(delta), zorder=depth, c=color, edgecolors=edge, alpha=1)
            ax1.plot(amps, deltas, zorder=2, alpha=0.5, color="k")
            ax2.plot(anchos, deltas, zorder=2, alpha=0.5, color="k")
    ax1.set_ylabel("$\delta$", fontsize=13)
    ax1.set_xlabel("$A_{\\textrm{eff}}$", fontsize=13)
    ax1.axes.tick_params(axis="both", labelsize=13)
    ax2.set_ylabel("$\delta$", fontsize=13)
    ax2.set_xlabel("$b_{\\textrm{eff}}$", fontsize=13)
    ax2.axes.tick_params(axis="both", labelsize=13)
    ax3.legend(fontsize=8)
    ax1.grid(alpha=0.2, zorder=1)
    ax2.grid(alpha=0.2, zorder=1)
    plt.tight_layout()
    plt.savefig("deltas_amplitudes.png", dpi=300)
    plt.close()

