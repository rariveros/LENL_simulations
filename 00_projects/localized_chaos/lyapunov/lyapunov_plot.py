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
    lyap_espectrums = []
    n = 0
    for directory in directories:
        directory = working_directory + "/" + directory
        files = os.listdir(directory)
        lyap_espectrum = np.loadtxt(directory + "/lyapunov_mean.txt", delimiter=',')
        D_ky = np.loadtxt(directory + "/D_ky.txt", delimiter=',')
        gamma = D_ky[0]
        lyap_espectrums.append(lyap_espectrum.tolist())
        lyap_espectrum = np.array(lyap_espectrum)
        Ies = np.arange(0, len(lyap_espectrum))
        lyap_espectrum = filtro_array(2, lyap_espectrum)
        plt.plot(Ies / len(lyap_espectrum), np.flip(np.sort(lyap_espectrum)), color=(1 - (n / 13), 0, n / 13), zorder=10, label="$\gamma_0 = " + f"{gamma:.3f}" + "$")
        n = n + 1
    plt.hlines(0, 0, 0.1, colors="k", zorder=0)
    plt.xlabel('$i/N$', size='20')
    plt.ylabel('$\\textrm{Lyapunov Spectrum}$', size='20')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlim(0, 0.06)
    plt.ylim([-0.11, 0.13])
    plt.grid(alpha=0.2)
    plt.legend()
    plt.tight_layout()
    plt.savefig('lyap_spectrum.png', dpi=300)
    plt.show()
    plt.close()