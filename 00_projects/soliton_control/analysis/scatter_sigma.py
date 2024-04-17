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
    sigmas = ["12.500", "15.000", "17.500"]
    gammas = ["0.180"]
    colors = (np.array([12.5, 15, 17.5]) - 12.5) / 5.01
    xd = 0
    fig = plt.figure()
    fig, ((ax1, ax2),(ax3, ax4)) = plt.subplots(nrows=2, ncols=2)
    if not os.path.exists(datafile):
        os.makedirs(datafile)
    for gamma in gammas:
        for i in range(len(sigmas)):
            sigma = sigmas[i]
            centers_i = []
            nus_i = []
            sigmas_i = []
            m = 0
            for directory in directories:
                m = m + 1
                if os.path.exists(working_directory + "/" + directory + '/sigma=' + sigma + '/gamma=' + gamma):
                    print("#############   " + directory + "   #############")
                    data = np.loadtxt(working_directory + "/" + directory + '/sigma=' + sigma + '/gamma=' + gamma + '/data.txt', delimiter=',')
                    params = np.loadtxt(working_directory + "/" + directory + '/sigma=' + sigma + '/gamma=' + gamma + '/parameters.txt', delimiter=',')
                    #[alpha, beta, gamma_0, mu_0, nu, sigma]
                    #[center, delta, ancho, amplitud]
                    nu = params[4]
                    nus_i.append(nu)
                    center = data[0]
                    centers_i.append(center)
                    sigmas_i.append(params[5])

                    print(m)
                    print(nu)
                    color = (colors[i], 0, 1 - colors[i])
                    depth = 10
                    edge = "black"
                    lab = "$\sigma_i=" + sigma + "$"
                    if len(data) == 5:
                        color = (1, 1, 1)
                        depth = 5
                        edge = "black"
                        lab = "xd"
                    if m == 5:
                        ax1.scatter(nu, np.abs(center), label=lab, zorder=depth, c=color, edgecolors=edge, alpha=1)
                    elif m == 20 and lab == "xd" and xd == 0:
                        ax1.scatter(nu, np.abs(center), label="$\\textrm{Oscillatory}$", zorder=depth, c=color, edgecolors=edge, alpha=1)
                        xd = 1
                    else:
                        ax1.scatter(nu, np.abs(center), zorder=depth, c=color, edgecolors=edge, alpha=1)
            ax1.plot(nus_i, np.abs(centers_i), zorder=3, color="k", alpha=1)
    ax1.set_xlabel("$\\nu$", fontsize=15)
    ax1.set_ylabel("$x_s$", fontsize=15)
    ax1.axes.tick_params(axis="both", labelsize=12)
    ax1.hlines(0, -0.17, -0.04, colors="k")
    ax1.set_xlim(-0.16, -0.05)
    ax1.legend()
    ax1.grid(alpha=0.2, zorder=0)
    plt.tight_layout()
    plt.savefig("parameter_space_sigma.png", dpi=300)
    plt.close()
