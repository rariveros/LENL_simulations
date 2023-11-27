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
    for directory_i in directories:
        directory = working_directory + "/" + directory_i
        files = os.listdir(directory)
        lyap_mean = np.loadtxt(directory + "/lyapunov_mean.txt", delimiter=',')
        D_ky = np.loadtxt(directory + "/D_ky.txt", delimiter=',')
        x = np.arange(len(lyap_mean)) / len(lyap_mean)
        D_kys.append(D_ky)
        plt.plot(x, np.flip(np.sort(lyap_mean)), label="$\\" + str(directory_i) + "$", linewidth=1)
    D_kys = np.array(D_kys)
    plt.title("$\\textrm{Lyapunov Spectrum}$", size='20')
    plt.legend(fontsize=8)
    plt.xlabel('$i/N$', size='20')
    plt.xticks([0.0, 0.1, 0.2], fontsize=15)
    plt.xlim([0, 0.1])
    plt.ylabel('$\lambda$', size='20')
    plt.yticks(fontsize=15)
    plt.ylim([-0.1, 0.1])
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig('lyap_espectrums.png', dpi=300)
    plt.close()
    #[gamma_0, sigma, D_ky]
    print(D_kys[:, 0])
    print(D_kys[:, 2])
    plt.plot(D_kys[:, 0], D_kys[:, 2], '-o')
    plt.xlabel('$\gamma_0$', size='20')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.ylabel('$D_{KY}$', size='20')
    plt.ylim([-2, 30])
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig('D_ky.png', dpi=300)
    plt.show()
    plt.close()