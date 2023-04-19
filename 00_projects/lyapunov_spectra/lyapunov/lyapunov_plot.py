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
    for directory_i in directories:
        directory = working_directory + "/" + directory_i
        files = os.listdir(directory)
        for file in files:
            if file[0] == "l":
                lyap_mean = np.loadtxt(directory + "/" + file, delimiter=',')
                x = np.arange(len(lyap_mean)) / len(lyap_mean)
            #elif file[0] == "x":
                #x_lyap = np.loadtxt(directory + "/" + file, delimiter=',')
        plt.plot(x, np.flip(np.sort(lyap_mean)), label="$\\" + str(directory_i) + "$", linewidth=1)
    plt.title("$\\textrm{Lyapunov Spectrum}$", size='20')
    plt.legend(fontsize=8)
    plt.xlabel('$i/N$', size='20')
    plt.xticks([0.0, 0.1, 0.2], fontsize=15)
    plt.xlim([0, 0.2])
    plt.ylabel('$\lambda$', size='20')
    plt.yticks(fontsize=15)
    #plt.ylim([-0.4, 1.5])
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig(initial_dir_data + '/lyap_espectrums_test_05.png', dpi=300)
    plt.close()