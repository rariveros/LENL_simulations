from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    frequencies = []
    disc = "aD:"
    initial_dir_data = str(disc) + 'mnustes_science/simulation_data/FD'
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    save_directory = working_directory
    directories = [name for name in os.listdir(working_directory) if os.path.isdir(os.path.join(working_directory, name))]
    MEAN_DEGREES = []
    N_QUIM = []
    N_QUIM_std = []
    for directory_01 in directories:
        dir_01 = working_directory + "/" + directory_01
        directories_01 = [name for name in os.listdir(dir_01) if os.path.isdir(os.path.join(dir_01, name))]
        n_quim = []
        mean_degree = []
        for directory_02 in directories_01:
            dir_02 = working_directory + "/" + directory_01 + "/" + directory_02
            output = np.loadtxt(dir_02 + '/output.txt', delimiter=',')
            n_quim.append(output[0])
            mean_degree.append(output[1])
        MEAN_DEGREES.append(mean_degree[0])
        N_QUIM.append(np.mean(n_quim))
        print(np.std(n_quim))
        N_QUIM_std.append(np.std(n_quim))
    MEAN_DEGREES = np.array(MEAN_DEGREES)
    N_QUIM = np.array(N_QUIM)
    N_QUIM_std = np.array(N_QUIM_std)

    fig, ax = plt.subplots(1, 1, figsize=(5, 3))
    ax.errorbar(MEAN_DEGREES, N_QUIM, yerr=N_QUIM_std, marker='o', ls='', ecolor="k", mec='black', color="k")
    ax.set_xlim(13, 33)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xticks([14, 16, 18, 20, 22, 24, 26, 28, 30, 32])
    ax.tick_params(axis="y", direction="in", labelsize=12, left=True, right=True, labelleft=True, labelright=False)
    ax.tick_params(axis="x", direction="in", labelsize=12, top=True, bottom=True, labeltop=False, labelbottom=True)
    ax.set_xlabel("$\\langle k \\rangle$", fontsize=20)
    ax.set_ylabel("$\\frac{N_c}{N}$", fontsize=20)
    ax.text(13.5, 0.92, "$\kappa = 1.5 \\times 10^{-2}$", fontsize=12)
    ax.text(13.5, 0.8, "$N = 501$", fontsize=12)
    plt.tight_layout()
    plt.savefig('mock_characterization.png', dpi=200)
