from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    disco = 'E:/'
    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FD'
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    files = [name for name in os.listdir(working_directory) if os.path.isfile(os.path.join(working_directory, name))]
    freqs = []
    ampds = []
    params = []
    sigmas = []
    power = []
    for file_i in files:
        if file_i[0:4] == "freq":
            freq_i = np.loadtxt(working_directory + '/' + file_i, delimiter=',')
            sigma_i = float(file_i.split("_")[-1].split(".")[0])
            freqs.append(freq_i)
            sigmas.append(sigma_i)
        elif file_i[0:4] == "ampd":
            ampd_i = np.loadtxt(working_directory + '/' + file_i, delimiter=',')
            ampds.append(ampd_i)
        elif file_i[0:4] == "para":
            param_i = np.loadtxt(working_directory + '/' + file_i, delimiter=',')
            params.append(param_i)
        elif file_i[0:4] == "powe":
            power_i = np.loadtxt(working_directory + '/' + file_i, delimiter=',')
            power.append(power_i)
    threshold = 0.002
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for i in range(len(freqs)):
        color = [(sigmas[i] - sigmas[0]) / (sigmas[-1] - sigmas[0]), 0, 1 - (sigmas[i] - sigmas[0]) / (sigmas[-1] - sigmas[0])]

        for j in range(len(params[i])):
            if power[i][j] ** 0.5 < threshold:
                color_j = "k"
                alpha_j = 1
                freqs[i][j] = 0
            else:
                color_j = color
                alpha_j = 1
            ax.scatter(sigmas[i], params[i][j], power[i][j] ** 0.5, c=color_j, alpha=alpha_j, edgecolors='black')
        ax.plot(sigmas[i] * np.ones(len(params[i])), params[i], power[i] ** 0.5, color=color, alpha=1)
    ax.xaxis.set_rotate_label(False)
    ax.yaxis.set_rotate_label(False)
    ax.zaxis.set_rotate_label(False)
    ax.set_xlabel('$\sigma_i$', fontsize=22, rotation=0)
    ax.set_ylabel('$\gamma_0$', fontsize=22, rotation=0)
    ax.set_zlabel('$P^{\\frac{1}{2}}$', fontsize=22, rotation=0)
    ax.axes.set_xlim3d(left=25, right=85)
    ax.axes.set_ylim3d(bottom=0.15, top=0.2)
    ax.axes.set_zlim3d(bottom=0)
    ax.axes.tick_params(axis="both", labelsize=12)
    ax.view_init(elev=21, azim=-145, roll=0)
    #plt.show()
    plt.tight_layout()
    plt.savefig(working_directory + '/power_3D.png', dpi=300)
    plt.close()