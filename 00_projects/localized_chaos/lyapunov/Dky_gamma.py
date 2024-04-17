import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *
from directories_lyap import *
from itertools import zip_longest

if __name__ == '__main__':
    disco = 'E:/'
    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FD' + main_directory + subdir_exponent
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    directories = os.listdir(working_directory)
    D_kys = []
    for directory_i in directories:
        directories_i = os.listdir(working_directory + "/" + directory_i)
        D_kys_i = []
        for directory in directories_i:
            directory = working_directory + "/" + directory_i + "/" + directory
            files = os.listdir(directory)
            D_ky = np.loadtxt(directory + "/D_ky.txt", delimiter=',')
            D_kys_i.append(D_ky.tolist())
        D_kys_i = np.array(D_kys_i)
        print(D_kys_i[0, :])
        sigma_i = D_kys_i[:, 1]
        plt.plot(D_kys_i[:, 0], D_kys_i[:, 2], 'o-', label="$\sigma_i = " + str(sigma_i[0]) + "$")
    plt.xlabel('$\gamma_0$', size='20')
    plt.xticks(fontsize=15)
    #plt.xlim(0.2, 0.525)
    plt.yticks(fontsize=15)
    plt.ylabel('$D_{KY}$', size='20')
    plt.ylim([-2, 30])
    plt.grid(alpha=0.2)
    plt.legend()
    plt.tight_layout()
    plt.savefig('D_ky_gamma.png', dpi=300)
    plt.show()
    plt.close()