import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    disco = 'E:/'
    initial_dir_data = str(disco) + 'Users/mnustes_science/PT_fluids/mnustes_science/simulation_data/simulation_dat'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    X = np.loadtxt(directory + '/X.txt', delimiter=',')
    T = np.loadtxt(directory + '/T.txt', delimiter=',')

    dx = X[1] - X[0]
    Nx = len(X)
    Z_complex = Z_r + 1j * Z_i
    arg = np.angle(Z_complex)
    Z_modulo = np.absolute(Z_complex)

    D = sparse_D_neumann(Nx, dx)
    dZ_r = Der(D, np.transpose(Z_r))
    dZ_i = Der(D, np.transpose(Z_i))
    dZ_modulo = dZ_r * dZ_i

    fig, ax = plt.subplots()

    # legend
    pcm = plt.pcolormesh(X, T, np.transpose(dZ_modulo), cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm)
    cbar.set_label('$A_R(x,t)$', rotation=0, size=25, labelpad=-27, y=1.11)

    # put the major ticks at the middle of each cell
    plt.xlabel('$x$', size='25')
    plt.xticks(fontsize=15)

    plt.ylabel('$t$', size='25')
    plt.yticks(fontsize=15)

    plt.grid(linestyle='--', alpha=0.2, color='k')
    cbar.ax.tick_params(labelsize=15)

    # labels
    plt.tight_layout()
    plt.show()
    #plt.savefig(directory + '/field_visual.png', dpi=300)
    #plt.close()