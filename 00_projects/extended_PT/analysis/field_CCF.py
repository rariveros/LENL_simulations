import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    disco = 'D:/'
    initial_dir_data = str(disco) + 'Users/mnustes_science/PT_fluids/mnustes_science/simulation_data'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    save_directory = directory + "//analysis"
    if not os.path.exists(save_directory):
        os.makedirs(save_directory)

    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    X = np.loadtxt(directory + '/X.txt', delimiter=',')
    T = np.loadtxt(directory + '/T.txt', delimiter=',')

    Nx = len(X)
    Nt = len(T)

    I_L = np.where(X == -10)[0][0]
    I_R = np.where(X == 10)[0][0]

    ti, tf = int(0.5 * Nt), Nt
    Z_r = Z_r[ti:tf, :]
    Z_i = Z_i[ti:tf, :]
    T = T[ti:tf] - T[ti]
    dt = T[1] - T[0]
    Nx = len(X)
    Nt = len(T)

    print(Nt)

    Z = Z_r + 1j * Z_i
    Z_conj = Z_r - 1j * Z_i
    Z_modulo = np.absolute(Z)

    plt.plot(T, Z_modulo[:, I_L], color="r", label="$|A(x_{I}, t)|$")
    plt.plot(T, Z_modulo[:, I_R], color="b", label="$|A(x_{R}, t)|$")
    plt.legend(fontsize="15", loc=1, framealpha=1)
    plt.grid(alpha=0.5)
    plt.xlabel('$t$', size='25')
    plt.xticks(fontsize=12)
    plt.xlim([T[0], T[-1]])
    plt.ylabel('$|A|$', size='25')
    plt.yticks(fontsize=12)
    plt.tight_layout()
    plt.savefig(save_directory + '/peak_oscillations.png', dpi=300)
    plt.close()

    CCF = []
    CCF_flipped = []

    for j in range(Nx):
        CCF_flipped_j = np.correlate(np.flip(Z, axis=1)[:, j], Z_conj[:, j] ,  "full")
        CCF_j = np.correlate(Z[:, j], Z_conj[:, j], "full")
        #CCF_flipped_j = np.correlate(np.flip(Z_modulo, axis=1)[:, j], Z_modulo[t_i_signal:t_f_signal, j],  "same")
        #CCF_j = np.correlate(Z_modulo[:, j], Z_modulo[t_i_signal:t_f_signal, j], "same")
        CCF_flipped.append(CCF_flipped_j)
        CCF.append(CCF_j)
    CCF = np.array(np.abs(CCF))
    CCF_flipped = np.array(np.abs(CCF_flipped))
    tau = np.arange(-T[-1], T[-1], dt) #np.arange(-T[-1], T[-1], dt)

    np.savetxt(save_directory + '/S+.txt', np.transpose(CCF), delimiter=',')
    np.savetxt(save_directory + '/S-.txt', np.transpose(CCF_flipped), delimiter=',')
    np.savetxt(save_directory + '/tau.txt', tau, delimiter=',')

    ################################

    plt.plot(tau, np.transpose(CCF)[:, I_L], label="$S_{+}(x_{R},\\tau)$", color=(216 / 255, 135 / 255, 25 / 255))
    plt.plot(tau, np.transpose(CCF_flipped)[:, I_L], label="$S_{-}(x_{R},\\tau)$", color=(49 / 255, 0 / 255, 129 / 255))
    plt.legend(fontsize="15", loc=1, framealpha=1)
    plt.grid(alpha=0.5)
    plt.xlabel('$\\tau$', size='25')
    plt.xticks(fontsize=12)
    plt.xlim([tau[0], tau[-1]])
    plt.ylabel('$|S_i(\\tau)|$', size='25')
    plt.yticks(fontsize=12)
    plt.tight_layout()
    plt.ylim(bottom=0)
    #plt.savefig(directory + '/peak_oscillations.png', dpi=300)
    plt.savefig(save_directory + '/peak_autocorrelation.png', dpi=300)
    plt.close()

    ################################

    # legend
    pcm = plt.pcolormesh(X, tau, np.transpose(CCF), cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm)
    cbar.set_label('$S_{+}(x, \\tau)$', rotation=0, size=25, labelpad=-27, y=1.11)

    # put the major ticks at the middle of each cell
    plt.xlabel('$x$', size='25')
    plt.xticks(fontsize=15)

    plt.ylabel('$\\tau$', size='25')
    plt.yticks(fontsize=15)
    #plt.ylim([-300, 300])

    plt.grid(linestyle='--', alpha=0.2, color='k')
    cbar.ax.tick_params(labelsize=15)

    # labels
    plt.tight_layout()
    plt.savefig(save_directory + '/S+.png', dpi=300)
    plt.close()

    ################################

    # legend
    pcm = plt.pcolormesh(X, tau, np.transpose(CCF_flipped), cmap="plasma", shading='auto')
    cbar = plt.colorbar(pcm)
    cbar.set_label('$S_{-}(x, \\tau)$', rotation=0, size=25, labelpad=-27, y=1.11)

    # put the major ticks at the middle of each cell
    plt.xlabel('$x$', size='25')
    plt.xticks(fontsize=15)

    plt.ylabel('$\\tau$', size='25')
    plt.yticks(fontsize=15)
    #plt.ylim([-300, 300])

    plt.grid(linestyle='--', alpha=0.2, color='k')
    cbar.ax.tick_params(labelsize=15)

    # labels
    plt.tight_layout()
    plt.savefig(save_directory + '/S-.png', dpi=300)
    plt.close()

    ################################

    # legend
    pcm = plt.pcolormesh(X, T, Z_modulo, cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm)
    cbar.set_label('$|A(x,t)|$', rotation=0, size=25, labelpad=-27, y=1.11)

    # put the major ticks at the middle of each cell
    plt.xlabel('$x$', size='25')
    plt.xticks(fontsize=15)

    plt.ylabel('$t$', size='25')
    plt.yticks(fontsize=15)
    #plt.ylim([-300, 300])

    plt.grid(linestyle='--', alpha=0.2, color='k')
    cbar.ax.tick_params(labelsize=15)

    # labels
    plt.tight_layout()
    plt.savefig(save_directory + '/module_+.png', dpi=300)
    plt.close()

    ################################

    # legend
    pcm = plt.pcolormesh(X, T, np.flip(Z_modulo, axis=1), cmap="plasma", shading='auto')
    cbar = plt.colorbar(pcm)
    cbar.set_label('$|A(-x,t)|$', rotation=0, size=25, labelpad=-27, y=1.11)

    # put the major ticks at the middle of each cell
    plt.xlabel('$x$', size='25')
    plt.xticks(fontsize=15)

    plt.ylabel('$t$', size='25')
    plt.yticks(fontsize=15)
    #plt.ylim([-300, 300])

    plt.grid(linestyle='--', alpha=0.2, color='k')
    cbar.ax.tick_params(labelsize=15)

    # labels
    plt.tight_layout()
    plt.savefig(save_directory + '/module_-.png', dpi=300)
    plt.close()

    ################################

    middle_tau = int(len(tau) / 2)
    plt.plot(X, np.transpose(CCF)[middle_tau, :], label="$S_{+}$", color=(216 / 255, 135 / 255, 25 / 255))
    plt.plot(X, np.transpose(CCF_flipped)[middle_tau, :], linestyle="--", label="$S_{-}$", color=(49 / 255, 0 / 255, 129 / 255))
    plt.legend(fontsize="15", loc=1, framealpha=1)
    plt.grid(alpha=0.5)
    plt.xlabel('$x$', size='25')
    plt.xticks(fontsize=12)
    plt.xlim([X[0], X[-1]])
    plt.ylabel('$S_{i}(x, \\tau_{\\textrm{max}})$', size='25')
    plt.yticks(fontsize=12)
    plt.tight_layout()
    plt.savefig(save_directory + '/S+_spatial.png', dpi=300)
    plt.close()