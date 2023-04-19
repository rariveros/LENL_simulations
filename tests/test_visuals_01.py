from functions import *
from back_process import *
from time_integrators import *
from directories_lyap import *

if __name__ == '__main__':
    disco = 'E:/'
    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FD' + main_directory
    root = tk.Tk()
    root.withdraw()
    directory_sim = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Carpeta de Simulación')
    Z_r = np.loadtxt(directory_sim + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory_sim + '/field_img.txt', delimiter=',')
    x_grid = np.loadtxt(directory_sim + '/X.txt', delimiter=',')
    t_grid = np.loadtxt(directory_sim + '/T.txt', delimiter=',')

    root = tk.Tk()
    root.withdraw()
    directory_ham = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Carpeta de Hamiltoniano')
    H_points = np.loadtxt(directory_ham + '/H_ext.txt', delimiter=',')
    T_points = np.loadtxt(directory_ham + '/T_ext.txt', delimiter=',')
    H_mean = np.loadtxt(directory_ham + '/H_mean.txt', delimiter=',')
    T_int = np.loadtxt(directory_ham + '/T_intervalos.txt', delimiter=',')
    T_ham = np.loadtxt(directory_ham + '/T.txt', delimiter=',')[0:-1]

    #################

    Z_complex = Z_r + 1j * Z_i
    Z_mod = np.absolute(Z_complex)

    fig = plt.figure()
    fig.set_figheight(6)
    fig.set_figwidth(9)
    ax1 = plt.subplot2grid(shape=(6, 9), loc=(0, 0), colspan=8, rowspan=2)
    ax2 = plt.subplot2grid(shape=(6, 9), loc=(2, 0), colspan=8, rowspan=4)
    cbar_ax = plt.subplot2grid(shape=(6, 18), loc=(2, 16), colspan=1, rowspan=4)
    axs = [ax1, ax2]

    # Some example data to display
    axs[0].plot(T_ham, H_mean, color="k")
    axs[0].scatter(T_points, H_points, c="r")
    axs[0].vlines(T_int, 0, 45, colors="r", alpha=0.5)
    axs[0].set_ylim([0, 45])
    axs[0].set_xlim([400, 1000])
    axs[0].tick_params(labelsize=10)
    #axs[0].set_xlabel('$t$', fontsize=15)
    axs[0].set_ylabel('$H\left[A, A^{*}\\right]$', fontsize=15)
    axs[0].grid(alpha=0.3, zorder=0)

    pcm = axs[1].pcolormesh(t_grid, x_grid, np.transpose(Z_mod), cmap=parula_map, shading='auto')
    axs[1].set_ylim([-25, 25])
    axs[1].set_xlim([400, 1000])
    axs[1].vlines(T_int, -25, 25, colors="r", alpha=1, zorder=1, linestyles="--")
    axs[1].tick_params(labelsize=10)
    axs[1].set_xlabel('$t$', fontsize=15)
    axs[1].set_ylabel('$x$', fontsize=15)
    axs[1].grid(alpha=0.3, zorder=0)
    # Make some room for the colorbar
    fig.subplots_adjust(left=0.1,
                        bottom=0.1,
                        right=0.9,
                        top=0.9,
                        wspace=1.2,
                        hspace=1.2)
    cbar = fig.colorbar(pcm, cax=cbar_ax)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label('$|A|$', rotation=0, size=15, labelpad=-27, y=1.1)

    plt.savefig(directory_sim + '/test.png', dpi=300)
    plt.close()