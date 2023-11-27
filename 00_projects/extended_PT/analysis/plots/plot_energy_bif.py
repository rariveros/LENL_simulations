from functions import *
from back_process import *
from time_integrators import *


if __name__ == '__main__':
    directory = "D:\mnustes_science\simulation_data\FD\PDNLS_extended_PT\dgamma\data_out"

    AR_up = np.loadtxt(directory + '/AR_up.txt', delimiter=',')
    AL_up = np.loadtxt(directory + '/AL_up.txt', delimiter=',')
    params_up = np.loadtxt(directory + '/params_up.txt', delimiter=',')
    AR_down = np.loadtxt(directory + '/AR_down.txt', delimiter=',')
    AL_down = np.loadtxt(directory + '/AL_down.txt', delimiter=',')
    params_down = np.loadtxt(directory + '/params_down.txt', delimiter=',')

    fig, ax = plt.subplots()
    ax.errorbar(params_up, AL_up[:, 0], AL_up[:, 1], marker='^', ls='', ecolor="k", mec='black',
                color=(216 / 255, 135 / 255, 25 / 255), label="$x_L$", alpha=0.8)
    ax.errorbar(params_down, AL_down[:, 0], AL_down[:, 1], marker='v', ls='', ecolor="k", mec='black',
                color=(49 / 255, 0 / 255, 129 / 255), label="$x_L$", alpha=0.8)
    plt.xlabel('$\gamma_0$', size='25')
    plt.ylabel('$A_{\\textrm{max}}$', size='25', rotation=90)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlim(0.17, 0.19)
    #plt.ylim(0, 2)
    #ax.yaxis.set_major_formatter(FormatStrFormatter('$%g \pi$'))
    #ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base=0.25))
    plt.grid(linestyle='--', alpha=0.5)
    #plt.legend(loc="upper right", fontsize=20)
    plt.tight_layout()
    plt.savefig(directory + '//bifurcation_L.png', dpi=300)
    plt.close()