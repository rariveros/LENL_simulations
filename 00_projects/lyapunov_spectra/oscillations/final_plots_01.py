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
            power_i = power[i, 1:7]
            params_i = params[i, 1:7]
        elif sigmas[i] == 75:
            power_i = power[i, 1:8]
            params_i = params[i, 1:8]
        elif sigmas[i] == 70:
            power_i = power[i, 1:8]
            params_i = params[i, 1:8]
        elif sigmas[i] == 65:
            power_i = power[i, 1:9]
            params_i = params[i, 1:9]
        elif sigmas[i] == 60:
            power_i = power[i, 1:10]
            params_i = params[i, 1:10]
        elif sigmas[i] == 55:
            power_i = power[i, 2:11]
            params_i = params[i, 2:11]
            ###########################
        elif sigmas[i] == 30:
            power_i = power[i, 10:]
            params_i = params[i, 10:]
        elif sigmas[i] == 35:
            power_i = power[i, 7:]
            params_i = params[i, 7:]
        elif sigmas[i] == 40:
            power_i = power[i, 6:]
            params_i = params[i, 6:]
        elif sigmas[i] == 45:
            power_i = power[i, 4:]
            params_i = params[i, 4:]
        elif sigmas[i] == 50:
            power_i = power[i, 3:]
            params_i = params[i, 3:]
            ###########################
        else:
            power_i = power[i]
            params_i = params[i]
        result = scipy.stats.linregress(params_i, power_i)
        regression_coef_i = - result.intercept / result.slope
        regression_coef.append(regression_coef_i)
        #plt.scatter(params_i, power_i ** 0.5)
    #plt.show()
    #plt.close()
    regression_coef = np.array(regression_coef)
    x = sigmas
    y = regression_coef
    spl = interpolate.splrep(x, y, s=1000)
    x2 = np.linspace(sigmas[0] - 5, sigmas[-1] + 5, 200)
    y2 = interpolate.splev(x2, spl)
    plt.plot(x2, y2, c="r", linestyle="--", linewidth=4)
    #plt.plot(x_int, y_int, c="k", zorder=10)
    plt.scatter(x, y, c="k", zorder=10, s=80)
    plt.grid(alpha=0.5, zorder=1)
    plt.xlim([28, 82])
    plt.ylim([0.15, 0.2])
    plt.xlabel("$\sigma_i$", fontsize=25)
    plt.ylabel("$\gamma_0$", fontsize=25)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tight_layout()
    plt.savefig('threshold.png', dpi=300)
    plt.close()