from functions import *
from back_process import *
from time_integrators import *
from directories_lyap import *

if __name__ == '__main__':
    disco = 'E:/'
    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FD' + main_directory
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    Z_r = np.loadtxt(directory + '/field_real_filt.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img_filt.txt', delimiter=',')
    X = np.loadtxt(directory + '/X.txt', delimiter=',')
    T = np.loadtxt(directory + '/T.txt', delimiter=',')

    Nx = len(X)
    Nt = len(T)
    dt = T[1] - T[0]
    dx = X[1] - X[0]
    L = X[-1] - X[0]

    Z_complex = Z_r + 1j * Z_i
    arg_light_1 = np.angle(Z_complex)
    arg_light_1 = (2 * np.pi + arg_light_1) * (arg_light_1 < 0) + arg_light_1 * (arg_light_1 > 0)
    phase_unwrapped = np.unwrap(arg_light_1)
    D = sparse_D(Nx, dx)
    k = Der(D, np.transpose(phase_unwrapped))
    k = np.transpose(k)
    pcm = plt.pcolormesh(X, T, np.abs(k), cmap=parula_map, shading='auto', vmin=0, vmax=4)
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$|k(x, t)|$', rotation=0, size=20, labelpad=-27, y=1.1)
    cbar.ax.tick_params(labelsize=15)
    plt.xlim([X[0], X[-1]])
    plt.xlabel('$x$', size='25')
    plt.ylabel('$t$', size='25')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.grid(linestyle='--', alpha=0.2, color='k')
    #plt.savefig(file + subfile + '/module_spacetime_01.png', dpi=300)
    plt.show()
    plt.close()