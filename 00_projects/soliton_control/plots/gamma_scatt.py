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
    gammas = ["0.100", "0.110", "0.120", "0.130", "0.140", "0.150", "0.160", "0.170", "0.180"]#, "0.190", "0.200", "0.210", "0.220", "0.240"]
    colors = (np.array([12.5, 15, 17.5]) - 12.5) / 5.01
    xd = 0
    fig = plt.figure()
    fig, ((ax1, ax2),(ax3, ax4)) = plt.subplots(nrows=2, ncols=2)
    if not os.path.exists(datafile):
        os.makedirs(datafile)
    centers_i = []
    centersR_i = []
    centersI_i = []
    gammas_i = []
    sigmas_i = []
    for directory in directories:
        print("#############   " + directory + "   #############")
        data = np.loadtxt(working_directory + "/" + directory + '/data.txt', delimiter=',')
        params = np.loadtxt(working_directory + "/" + directory + '/parameters.txt', delimiter=',')
        #[alpha, beta, gamma_0, mu_0, nu, sigma]
        #[center, delta, ancho, amplitud]
        gamma_0 = params[2]
        gammas_i.append(gamma_0)
        center = data[0]
        centers_i.append(center)
        centerR = data[1]
        centersR_i.append(centerR)
        centerI = data[2]
        centersI_i.append(centerI)
        sigmas_i.append(params[5])
        print(gamma_0)

        depth = 10
        edge = "black"
        ax1.scatter(gamma_0, np.abs(center), zorder=depth, c="k", edgecolors=edge, alpha=1)
        ax1.scatter(gamma_0, np.abs(centerR), zorder=depth, c="b", edgecolors=edge, alpha=1)
        ax1.scatter(gamma_0, np.abs(centerI), zorder=depth, c="r", edgecolors=edge, alpha=1)
    ax1.plot(gammas_i, np.abs(centers_i), zorder=3, color="k", alpha=1)
    ax1.plot(gammas_i, np.abs(centersR_i), zorder=3, color="k", alpha=1)
    ax1.plot(gammas_i, np.abs(centersI_i), zorder=3, color="k", alpha=1)
    ax1.set_xlabel("$\gamma$", fontsize=15)
    ax1.set_ylabel("$x_s$", fontsize=15)
    ax1.axes.tick_params(axis="both", labelsize=12)
    #ax1.hlines(0, -0.17, -0.04, colors="k")
    #ax1.set_xlim(-0.16, -0.05)
    ax1.grid(alpha=0.2, zorder=0)
    plt.tight_layout()
    plt.savefig("parameter_space_sigma.png", dpi=300)
    plt.close()
