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
    sigmas = ["15.000"]
    gammas = ["0.100", "0.110", "0.120", "0.130", "0.140", "0.150", "0.160", "0.170", "0.180", "0.190", "0.200", "0.210", "0.220", "0.240"]
    colors = (np.array([0.16, 0.18, 0.20]) - 0.16) / 0.04
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    if not os.path.exists(datafile):
        os.makedirs(datafile)
    for gamma in gammas:
        for i in range(len(sigmas)):
            sigma = sigmas[i]
            centers_i = []
            nus_i = []
            sigmas_i = []
            for directory in directories:
                if os.path.exists(working_directory + "/" + directory + '/sigma=' + sigma + '/gamma=' + gamma):
                    print("#############   " + directory + "   #############")
                    data = np.loadtxt(working_directory + "/" + directory + '/sigma=' + sigma + '/gamma=' + gamma + '/data.txt', delimiter=',')
                    params = np.loadtxt(working_directory + "/" + directory + '/sigma=' + sigma + '/gamma=' + gamma + '/parameters.txt', delimiter=',')
                    #[alpha, beta, gamma_0, mu_0, nu, sigma]
                    #[center, delta, ancho, amplitud]
                    nu = params[4]
                    center = data[0]
                    centers_i.append(center)
                    nus_i.append(nu)
                    sigmas_i.append(params[5])
                    color = (colors[i], 0, 1 - colors[i])
                    depth = 10
                    edge = "black"
                    if len(data) == 5:
                        color = (1, 1, 1)
                        depth = 1
                        edge = "white"
                    ax.scatter(nu, params[5], np.abs(center), label="$\sigma_i=" + sigma + "$", zorder=depth, c=color, s=40, edgecolors=edge, alpha=1)
            np.savetxt(working_directory + "/" + directory + '/centers.txt', centers_i, delimiter=',')
            np.savetxt(working_directory + "/" + directory + '/nus.txt', nus_i, delimiter=',')
            ax.plot(nus_i, sigmas_i, np.abs(centers_i), label="$\sigma_i=" + sigma + "$", zorder=3, color="k", alpha=1)

    plt.xlabel("$\\nu$", fontsize=20)
    plt.ylabel("$\sigma_i$", fontsize=20)
    ax.set_zlabel("$\Delta x$", fontsize=20)
    ax.view_init(elev=29, azim=53, roll=0)
    ax.axes.tick_params(axis="both", labelsize=15)
    #plt.legend()
    plt.grid(alpha=0.2, zorder=0)
    plt.tight_layout()
    plt.savefig("xs_gamma.png", dpi=300)
    plt.close()

