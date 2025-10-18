from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    working_directory = r"D:\mnustes_science\simulation_data\FD\coherent\phase_space\Delta=0.0000\gamma=0.1000\Omega=0.2000\g=0.0010"
    directories = [name for name in os.listdir(working_directory) if os.path.isdir(os.path.join(working_directory, name))]
    dir_0 = working_directory + "/" + directories[0]
    dir_1 = working_directory + "/" + directories[1]
    dir_2 = working_directory + "/" + directories[2]

    H01 = np.loadtxt(dir_0 + "/histogram_counts1.txt")
    H02 = np.loadtxt(dir_0 + "/histogram_counts2.txt")
    H11 = np.loadtxt(dir_1 + "/histogram_counts1.txt")
    H12 = np.loadtxt(dir_1 + "/histogram_counts2.txt")
    H21 = np.loadtxt(dir_2 + "/histogram_counts1.txt")
    H22 = np.loadtxt(dir_2 + "/histogram_counts2.txt")

    x01 = np.loadtxt(dir_0 + "/xedges1.txt")
    x02 = np.loadtxt(dir_0 + "/xedges2.txt")
    x11 = np.loadtxt(dir_1 + "/xedges1.txt")
    x12 = np.loadtxt(dir_1 + "/xedges2.txt")
    x21 = np.loadtxt(dir_2 + "/xedges1.txt")
    x22 = np.loadtxt(dir_2 + "/xedges2.txt")

    y01 = np.loadtxt(dir_0 + "/yedges1.txt")
    y02 = np.loadtxt(dir_0 + "/yedges2.txt")
    y11 = np.loadtxt(dir_1 + "/yedges1.txt")
    y12 = np.loadtxt(dir_1 + "/yedges2.txt")
    y21 = np.loadtxt(dir_2 + "/yedges1.txt")
    y22 = np.loadtxt(dir_2 + "/yedges2.txt")

    X01, Y01 = np.meshgrid(x01, y01)
    X02, Y02 = np.meshgrid(x02, y02)
    X11, Y11 = np.meshgrid(x11, y11)
    X12, Y12 = np.meshgrid(x12, y12)
    X21, Y21 = np.meshgrid(x21, y21)
    X22, Y22 = np.meshgrid(x22, y22)

    fig, ((ax01, ax11, ax21), (ax02, ax12, ax22)) = plt.subplots(2, 3, figsize=(4, 3))
    ax01.pcolormesh(X01, Y01, H01.T, cmap="afmhot")
    ax01.set_xlabel(r"$\kappa = 0.40$", fontsize=15, labelpad=15)
    ax01.set_ylabel("$\\textrm{Im }(\\alpha_j)$", fontsize=15)
    #ax01.set_xticks([-1, 0, 1])
    #ax01.set_yticks([-2, 0, 2])
    ax01.tick_params(axis="both", direction="in", labelsize=11, top=True, bottom=True, labeltop=True, labelbottom=False)
    ax01.xaxis.set_label_position("top")

    ax11.pcolormesh(X11, Y11, H11.T, cmap="afmhot")
    ax11.set_xlabel(r"$\kappa = 0.51$", fontsize=22, labelpad=15)
    #ax11.set_xticks([-1, 0, 1])
    #ax11.set_yticks([-2, 0, 2])
    ax11.tick_params(axis="x", direction="in", labelsize=11, top=True, bottom=True, labeltop=True, labelbottom=False)
    ax11.tick_params(axis="y", direction="in", labelsize=11, left=True, right=True, labelleft=False, labelright=False)
    ax11.xaxis.set_label_position("top")

    ax21.pcolormesh(X21, Y21, H21.T, cmap="afmhot")
    ax21.set_xlabel(r"$\kappa = 0.60$", fontsize=15, labelpad=15)
    #ax21.set_xticks([-1, 0, 1])
    #ax21.set_yticks([-2, 0, 2])
    ax21.tick_params(axis="both", direction="in", labelsize=11, top=True, bottom=True, labeltop=True, labelbottom=False)
    ax21.tick_params(axis="y", direction="in", labelsize=11, left=True, right=True, labelleft=False, labelright=False)
    ax21.xaxis.set_label_position("top")

    ### Fila 2 ###

    ax02.pcolormesh(X02, Y02, H02.T, cmap="afmhot")
    ax02.set_xlabel(r"$\kappa = 0.40$", fontsize=18)
    ax02.set_ylabel("$\\textrm{Im }(\\alpha_j)$", fontsize=15)
    #ax02.set_xticks([-2, 0, 2])
    #ax02.set_yticks([-1, 0, 1])
    ax02.tick_params(axis="both", direction="in", labelsize=11, top=True, bottom=True, labeltop=False, labelbottom=True)

    ax12.pcolormesh(X12, Y12, H12.T, cmap="afmhot")
    ax12.set_xlabel(r"$\kappa = 0.51$", fontsize=18)
    #ax12.set_xticks([-2, 0, 2])
    #ax12.set_yticks([-1, 0, 1])
    ax12.tick_params(axis="x", direction="in", labelsize=11, top=True, bottom=True, labeltop=False, labelbottom=True)
    ax12.tick_params(axis="y", direction="in", labelsize=11, left=True, right=True, labelleft=False, labelright=False)

    ax22.pcolormesh(X22, Y22, H22.T, cmap="afmhot")
    #ax22.set_xlabel("$\\textrm{Re }(\\alpha_j)$", fontsize=
    ax22.set_xlabel(r"$\kappa = 0.60$", fontsize=18)
    #ax22.set_xticks([-2, 0, 2])
    #ax22.set_yticks([-1, 0, 1])
    ax22.tick_params(axis="both", direction="in", labelsize=11, top=True, bottom=True, labeltop=False, labelbottom=True)
    ax22.tick_params(axis="y", direction="in", labelsize=11, left=True, right=True, labelleft=False, labelright=False)

    fig.subplots_adjust(wspace=0.05, hspace=0.05, left=0.15, bottom=0.2, right=0.99) #left=0.1, right=0.9, bottom=0.25, top=0.8)

    plt.savefig("phase_space.png", dpi=300)
    plt.close()
