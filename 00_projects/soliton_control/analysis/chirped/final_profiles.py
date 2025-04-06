import matplotlib.pyplot as plt
from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    directory = "C:/mnustes_science/simulation_data/FD/chirped_soliton/alpha=6.524/beta=1.000/mu=0.075/nu=-0.130/sigma=15.000/gamma=0.180"
    mods = np.loadtxt(directory + '/final_profiles.txt', delimiter=',')
    x_grid = np.loadtxt(directory + '/X.txt', delimiter=',')
    Cs = np.loadtxt(directory + '/Cs.txt', delimiter=',')

    pcm = plt.pcolormesh(x_grid, Cs, mods, cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$|A|$', rotation=0, size=20, labelpad=-27, y=1.1)
    cbar.ax.tick_params(labelsize=15)
    #plt.xlim(-40, 40)
    plt.xlabel('$x$', size='25')
    plt.ylabel('$c$', size='25')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.grid(linestyle='--', alpha=0.2, color='k')
    plt.tight_layout()
    plt.show()