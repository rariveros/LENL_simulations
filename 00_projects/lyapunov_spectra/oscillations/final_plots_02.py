from functions import *
from back_process import *
from time_integrators import *
import scipy
from scipy import interpolate

if __name__ == '__main__':
    disco = 'E:/'
    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FDa'
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
    power = np.array(power)
    params = np.array(params)
    sigmas = np.array(sigmas)
    gammas = params[0]
    chaos_or_osc = np.ones((len(gammas), len(sigmas)))
    print(chaos_or_osc)
    regression_coef = []
    for i in range(len(sigmas)):
        if sigmas[i] == 80:
            params_osc_i = params[i, 1:8]
            params_chaos_i = params[i, 8:]
            params_stat_i = np.array(params[i, 0:1])
        elif sigmas[i] == 75:
            params_osc_i =params[i, 1:9]
            params_chaos_i = params[i, 9:]
            params_stat_i = np.array(params[i, 0:1])
        elif sigmas[i] == 70:
            params_osc_i = params[i, 1:9]
            params_chaos_i = params[i, 9:]
            params_stat_i = np.array(params[i, 0:1])
        elif sigmas[i] == 65:
            params_osc_i = params[i, 2:10]
            params_chaos_i = params[i, 10:]
            params_stat_i = np.array(params[i, 0:2])
        elif sigmas[i] == 60:
            params_osc_i = params[i, 2:11]
            params_chaos_i = params[i, 11:]
            params_stat_i = np.array(params[i, 0:2])
        elif sigmas[i] == 55:
            params_osc_i = params[i, 3:12]
            params_chaos_i = params[i, 12:]
            params_stat_i = params[i, 0:3]
            ###########################
        elif sigmas[i] == 30:
            params_osc_i = params[i, 10:]
            params_stat_i = params[i, 0:10]
        elif sigmas[i] == 35:
            params_osc_i = params[i, 8:]
            params_stat_i = params[i, 0:8]
        elif sigmas[i] == 40:
            params_osc_i = params[i, 6:]
            params_stat_i = params[i, 0:6]
        elif sigmas[i] == 45:
            params_osc_i = params[i, 5:]
            params_stat_i = params[i, 0:5]
        elif sigmas[i] == 50:
            params_osc_i = params[i, 4:]
            params_stat_i = params[i, 0:4]
            ###########################
        else:
            params_stat_i = params[i]
        plt.scatter(sigmas[i] * np.ones(len(params_osc_i)), params_osc_i, facecolors='r', edgecolors='k', zorder=10)
        plt.scatter(sigmas[i] * np.ones(len(params_stat_i)), params_stat_i, facecolors='none', edgecolors='k', zorder=10)
        if sigmas[i] > 52:
            plt.scatter(sigmas[i] * np.ones(len(params_chaos_i)), params_chaos_i, c="k", marker="x", zorder=10)

    plt.grid(alpha=0.5, zorder=0)
    plt.xlabel("$\sigma_i$", fontsize=25)
    plt.ylabel("$\gamma_0$", fontsize=25)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tight_layout()
    plt.savefig('map.png', dpi=300)
    plt.close()