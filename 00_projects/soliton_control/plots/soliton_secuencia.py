from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    main_directory = "D:/mnustes_science/simulation_data/FD/soliton_control/alpha=6.524/beta=1.000/mu=0.075"
    directory_01 = "/nu=-0.150/sigma=15.000/gamma=0.180"
    directory_02 = "/nu=-0.200/sigma=15.000/gamma=0.180"
    directory_03 = "/nu=-0.170/sigma=15.000/gamma=0.180"

    params_01 = np.loadtxt(main_directory + "/" + directory_01 + '/parameters.txt', delimiter=',')
    ZR_01 = np.loadtxt(main_directory + "/" + directory_01 + '/field_real.txt', delimiter=',')
    ZI_01 = np.loadtxt(main_directory + "/" + directory_01 + '/field_img.txt', delimiter=',')
    Z_01 = ZR_01 + 1j * ZI_01
    Z_mod_01 = np.abs(Z_01)
    X_01 = np.loadtxt(main_directory + "/" + directory_01 + '/X.txt', delimiter=',')
    T_01 = np.loadtxt(main_directory + "/" + directory_01 + '/T.txt', delimiter=',')

    params_02 = np.loadtxt(main_directory + "/" + directory_02 + '/parameters.txt', delimiter=',')
    ZR_02 = np.loadtxt(main_directory + "/" + directory_02 + '/field_real.txt', delimiter=',')
    ZI_02 = np.loadtxt(main_directory + "/" + directory_02 + '/field_img.txt', delimiter=',')
    Z_02 = ZR_02 + 1j * ZI_02
    Z_mod_02 = np.abs(Z_02)
    X_02 = np.loadtxt(main_directory + "/" + directory_02 + '/X.txt', delimiter=',')
    T_02 = np.loadtxt(main_directory + "/" + directory_02 + '/T.txt', delimiter=',')

    params_03 = np.loadtxt(main_directory + "/" + directory_03 + '/parameters.txt', delimiter=',')
    ZR_03 = np.loadtxt(main_directory + "/" + directory_03 + '/field_real.txt', delimiter=',')
    ZI_03 = np.loadtxt(main_directory + "/" + directory_03 + '/field_img.txt', delimiter=',')
    Z_03 = ZR_03 + 1j * ZI_03
    Z_mod_03 = np.abs(Z_03)
    X_03 = np.loadtxt(main_directory + "/" + directory_03 + '/X.txt', delimiter=',')
    T_03 = np.loadtxt(main_directory + "/" + directory_03 + '/T.txt', delimiter=',')

    Z = np.vstack((Z_mod_01, Z_mod_02, Z_mod_03))
    dt = T_01[1] - T_01[0]
    T = np.arange(0, len(Z)) * dt
    plt.plot(T)
    plt.show()
    x_grid = X_01

    pcm = plt.pcolormesh(x_grid, T, Z, cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$|A|$', rotation=0, size=20, labelpad=-27, y=1.1)
    cbar.ax.tick_params(labelsize=15)
    plt.xlim([x_grid[0], x_grid[-1]])
    plt.xlabel('$x$', size='25')
    plt.ylabel('$t$', size='25')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    # plt.xlim([-50, 50])
    plt.grid(linestyle='--', alpha=0.2, color='k')
    plt.tight_layout()
    plt.show()
    #plt.savefig(file + subfile + '/module_spacetime.png', dpi=100)
    #plt.close()