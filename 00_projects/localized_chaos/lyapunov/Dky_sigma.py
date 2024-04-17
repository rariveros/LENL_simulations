import matplotlib.pyplot as plt
from functions import *
from back_process import *
from time_integrators import *
from directories_lyap import *
from itertools import zip_longest
from scipy import stats

if __name__ == '__main__':
    disco = 'E:/'
    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FD' + main_directory + subdir_exponent
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    directories = os.listdir(working_directory)
    D_kys = []
    sigmas = []
    gammas = []
    for directory_i in directories:
        directories_i = os.listdir(working_directory + "/" + directory_i)
        D_kys_i = []
        for directory in directories_i:
            directory = working_directory + "/" + directory_i + "/" + directory
            files = os.listdir(directory)
            D_ky = np.loadtxt(directory + "/D_ky.txt", delimiter=',')
            D_kys_i.append(D_ky.tolist())
        D_kys_i = np.array(D_kys_i)
        sigma_i = D_kys_i[:, 1]
        sigmas.append(D_kys_i[-5:, 1])
        D_kys.append(D_kys_i[-5:, 2])
        gammas.append(D_kys_i[-5:, 0])
        #plt.plot(D_kys_i[-:, 1], D_kys_i[-4:, 2], 'o', label="$\sigma_i = " + str(sigma_i[0]) + "$")
    sigmas = np.array(sigmas)
    D_kys = np.array(D_kys)
    deltas = []
    gamma = []
    n = 8
    for i in range(len(sigmas[0, :])):
        sigma_i = np.sort(sigmas[:, i])
        D_ky_i = np.sort(D_kys[:, i])
        slope, intercept, r, p, se = stats.linregress(sigma_i, D_ky_i)
        print(sigmas[:, i])
        print(D_kys[:, i])
        deltas.append(1/slope)
        gamma.append(gammas[0][i])
        plt.plot(sigmas[:, i], slope * sigmas[:, i] + intercept, linestyle="--", color=(1 - (n / 13), 0, n / 13), label="$\gamma_0 = " + f"{gammas[0][i]:.3f}" + "$")
        plt.scatter(sigmas[:, i], D_kys[:, i], color=(1 - (n / 13), 0, n / 13),)
        n = n + 1
    plt.xlabel('$\sigma_i$', size='20')
    plt.xticks(fontsize=15)
    #plt.xlim(0.2, 0.525)
    plt.yticks(fontsize=15)
    plt.ylabel('$D_{KY}$', size='20')
    plt.ylim([-2, 30])
    plt.grid(alpha=0.2)
    plt.legend()
    plt.tight_layout()
    plt.savefig('D_ky_sigma.png', dpi=300)
    #plt.show()
    plt.close()

    plt.scatter(gamma, deltas, s=30, c="k")
    plt.xlabel('$\gamma_0$', size='20')
    plt.xticks(fontsize=15)
    #plt.xlim(0.2, 0.525)
    plt.yticks(fontsize=15)
    plt.ylabel('$\zeta_{\delta}$', size='20')
    plt.ylim([0, 0.5])
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig('delta_Dky.png', dpi=300)
    plt.show()
    plt.close()

