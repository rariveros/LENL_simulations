import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    disco = 'E:/'
    initial_dir_data = str(disco) + "mnustes_science/simulation_data/FD/PDNLS_chaos/branching_dynamicsa"
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    Z_r_fft = np.loadtxt(working_directory + '/field_real_FFT.txt', delimiter=',')
    X = np.loadtxt(working_directory + '/K.txt', delimiter=',')
    T = np.loadtxt(working_directory + '/T.txt', delimiter=',')
    T = T
    Nx = len(X)
    Nt = len(T)
    #Z_filtered = filtro_superficie(Z_r_fft, 150, "X")

    #fig, ax = plt.subplots()

    # legend
    #pcm = plt.pcolormesh(X, T, Z_filtered, cmap="jet", shading='auto')
    #cbar = plt.colorbar(pcm)
    #cbar.set_label('$\hat{A}_R(t, k)$', rotation=0, size=20, labelpad=-27, y=1.11)

    # put the major ticks at the middle of each cell
    #plt.xlim([0, 0.6])
    #plt.xlabel('$k$', size='20')
    #plt.xticks(fontsize=15)

    #plt.ylabel('$T$', size='20')
    #plt.yticks(fontsize=15)
    # plt.ylim([0, 0.5])

    #plt.grid(linestyle='--', alpha=0.2, color='k')
    #cbar.ax.tick_params(labelsize=15)

    # labels
    #plt.tight_layout()
    #plt.show()
    #plt.close()
    K_max = []
    for i in range(Nt):
        K_max_i = X[np.argmax(Z_r_fft[i, :])]
        K_max.append(K_max_i)
    K_max = np.array(K_max)
    plt.scatter(T, K_max, s=1, c="r")
    plt.plot(T, K_max, color="k")
    plt.xlim([T[0], T[-1]])
    plt.grid(alpha=0.2, color="k")
    plt.show()
    plt.close()
