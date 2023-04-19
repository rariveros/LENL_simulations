import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *
from directories_lyap import *

if __name__ == '__main__':
    disco = 'E:/'
    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FD' + main_directory + subdir_exponent
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    directories = os.listdir(working_directory)
    D_kys = []
    gammas = []
    for directory_i in directories:
        directory = working_directory + "/" + directory_i
        files = os.listdir(directory)
        for file in files:
            if file[0] == "l":
                lyap_mean = np.loadtxt(directory + "/" + file, delimiter=',')
                x = np.arange(len(lyap_mean)) / len(lyap_mean)
                sum = 0
                if lyap_mean[0] <= 0:
                    D_ky = 0
                else:
                    for i in range(len(x)):
                        if sum >= 0:
                            sum = sum + lyap_mean[i]
                            first_neg = lyap_mean[i + 1]
                            p = i
                    D_ky = p + sum / first_neg
                D_kys.append(D_ky)
                gamma = directory_i.split("=")
                gamma = float(gamma[-1])
                gammas.append(gamma)
    D_kys = np.array(D_kys)
    gammas = np.array(gammas)
    save_directory = initial_dir_data + '/test_05'
    if not os.path.exists(save_directory):
        os.makedirs(save_directory)
    np.savetxt(save_directory + '/KY_dim_sigma.txt', D_kys, delimiter=',')
    np.savetxt(save_directory + '/sigmas.txt', gammas, delimiter=',')

    plt.scatter(gammas, D_kys, color="k")
    plt.xlabel('$\gamma_0$', size='20')
    plt.ylabel('$D_{KY}$', size='20')
    plt.xlim([5, 1.1 * np.amax(gammas)])
    plt.ylim([0, 1.1 * np.amax(D_kys)])
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig(initial_dir_data + '/kaplan-york_test_05.png', dpi=300)
    plt.close()