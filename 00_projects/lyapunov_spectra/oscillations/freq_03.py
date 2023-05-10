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
            sigma_i = float(file_i[-6:-4])
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
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for i in range(len(freqs)):
        ax.scatter(sigmas[i], params[i], power[i] ** 0.5)
        ax.plot(sigmas[i] * np.ones(len(params[i])), params[i], power[i] ** 0.5)
    ax.set_xlabel('$\sigma_i$', fontsize=20)
    ax.set_ylabel('$\gamma_0$', fontsize=20)
    ax.set_zlabel('$\\textrm{Average Power}$', fontsize=20)
    plt.show()