from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    NUS = np.loadtxt('C:/Users/Rafa/PycharmProjects/LENL_simulations/00_projects/rabi_windows/figures/NUS.txt', delimiter=',')
    FREQ = np.loadtxt('C:/Users/Rafa/PycharmProjects/LENL_simulations/00_projects/rabi_windows/figures/FREQS.txt', delimiter=',')
    distances = np.loadtxt('C:/Users/Rafa/PycharmProjects/LENL_simulations/00_projects/rabi_windows/figures/DIST.txt', delimiter=',')

    FREQ = np.array(FREQ)
    FREQ[np.isnan(FREQ)] = 0
    FREQ = filtro_superficie(FREQ, 4, "YX")

    pcm = plt.pcolormesh(distances, np.array(NUS), np.array(FREQ), cmap="jet", shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    plt.xlabel('$d$', size='20')
    plt.ylabel('$\\nu$', size='20')
    cbar.set_label('$\Omega \\times 10^{-3}$', rotation=0, size=15, labelpad=-27, y=1.1)
    plt.grid(linestyle='--', alpha=0.5)
    plt.savefig("freqs_var_nu.png", dpi=250)
    plt.close()
