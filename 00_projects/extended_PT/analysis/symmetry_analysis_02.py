import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *


if __name__ == '__main__':
    disco = 'D:/'
    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FDa'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    print("#############   " + directory + "   #############")
    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    X = np.loadtxt(directory + '/X.txt', delimiter=',')
    T = np.loadtxt(directory + '/T.txt', delimiter=',')
    params = np.loadtxt(directory + '/parameters.txt', delimiter=',')
    center_L = - 0.5 * params[3]
    center_R = 0.5 * params[3]

    Nx = len(X)
    Nt = len(T)

    I_L = np.where(X == center_L)[0][0]
    I_R = np.where(X == center_R)[0][0]

    ti, tf = int(0.4 * Nt), Nt
    Z_r = Z_r[ti:tf, :]
    Z_i = Z_i[ti:tf, :]
    T = T[ti:tf] - T[ti]
    dt = T[1] - T[0]

    Nx = len(X)
    Nt = len(T)

    Z = Z_r + 1j * Z_i
    Z_conj = Z_r - 1j * Z_i
    Z_modulo = np.absolute(Z)

    Z_mod_L_points_max, tau_R_points_max, I_R_max = max_finder(Z_modulo[:, I_L], T, Nt, dt)
    Z_mod_R_points_max, tau_L_points_max, I_L_max = max_finder(Z_modulo[:, I_R], T, Nt, dt)
    Z_sim_r = np.flip(Z_r, axis=1)
    Z_sim_i = np.flip(Z_i, axis=1)

    Z_sim = Z_sim_r + 1j * Z_sim_i
    Z_sim = 1j * Z_sim

    Z_sim_r = np.real(Z_sim)
    Z_sim_i = np.imag(Z_sim)

    plt.subplot(2, 1, 1)
    #plt.title("$A_R(x,t)$")
    plt.plot(X, Z_r[I_R_max[0], :], label="$\\textrm{Re}\left\{A(x, t)\\right\}$", c="b")
    plt.plot(X, Z_i[I_R_max[0], :], label="$\\textrm{Im}\left\{A(x, t)\\right\}$", c="r")
    plt.legend(loc="upper right", fontsize=12)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlim([X[0], X[-1]])
    plt.xlabel("$x$", fontsize=30)
    plt.grid(alpha=0.3)

    plt.subplot(2, 1, 2)
    #plt.title("$A_I(x,t)$")
    plt.plot(X, Z_sim_r[I_L_max[0], :], label="$\\textrm{Re}\left\{\mathcal{P}e^{i\pi/2}A(x, t + \\tau/2)\\right\}$", c="b")
    plt.plot(X, Z_sim_i[I_L_max[0], :], label="$\\textrm{Im}\left\{\mathcal{P}e^{i\pi/2}A(x, t + \\tau/2)\\right\}$", c="r")
    plt.legend(loc="upper right", fontsize=12)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlim([X[0], X[-1]])
    plt.xlabel("$x$", fontsize=30)
    plt.grid(alpha=0.3)

    plt.tight_layout()
    #plt.show()
    plt.savefig(directory + '/diff_1d.png', dpi=300)
    plt.close()

    # legend
    pcm = plt.pcolormesh(X, T[I_R_max[0]:Nt - (I_L_max[0] + 3)], Z_r[I_R_max[0]:Nt - (I_L_max[0] + 3), :], cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm)
    cbar.set_label('$A_R(x, t)$', rotation=0, size=25, labelpad=-27, y=1.11)

    # put the major ticks at the middle of each cell
    plt.xlim([X[0], X[-1]])
    plt.xlabel('$x$', size='25')
    plt.xticks(fontsize=15)

    plt.ylim([0, tf - ti])
    plt.ylabel('$t$', size='25')
    plt.yticks(fontsize=15)

    plt.grid(linestyle='--', alpha=0.2, color='k')
    cbar.ax.tick_params(labelsize=15)

    # labels
    plt.tight_layout()
    plt.savefig(directory + '/field_r.png', dpi=300)
    plt.close()

    # legend
    pcm = plt.pcolormesh(X, T[I_R_max[0]:Nt - (I_L_max[0] + 3)], Z_i[I_R_max[0]:Nt - (I_L_max[0] +3), :], cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm)
    cbar.set_label('$A_I(x, t)$', rotation=0, size=25, labelpad=-27, y=1.11)

    # put the major ticks at the middle of each cell
    plt.xlim([X[0], X[-1]])
    plt.xlabel('$x$', size='25')
    plt.xticks(fontsize=15)

    plt.ylim([0, tf - ti])
    plt.ylabel('$t$', size='25')
    plt.yticks(fontsize=15)

    plt.grid(linestyle='--', alpha=0.2, color='k')
    cbar.ax.tick_params(labelsize=15)

    # labels
    plt.tight_layout()
    plt.savefig(directory + '/field_i.png', dpi=300)
    plt.close()

    # --------------- #

    # legend
    pcm = plt.pcolormesh(X, T[I_R_max[0]:Nt - (I_L_max[0] + 3)], Z_sim_r[(I_L_max[0] + 3):Nt - I_R_max[0], :],
                         cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm)
    cbar.set_label('$A\'_R(x,t\')$', rotation=0, size=25, labelpad=-27, y=1.11)

    # put the major ticks at the middle of each cell
    plt.xlim([X[0], X[-1]])
    plt.xlabel('$x$', size='25')
    plt.xticks(fontsize=15)

    plt.ylim([0, tf - ti])
    plt.ylabel('$t\'$', size='25')
    plt.yticks(fontsize=15)

    plt.grid(linestyle='--', alpha=0.2, color='k')
    cbar.ax.tick_params(labelsize=15)

    # labels
    plt.tight_layout()
    plt.savefig(directory + '/Ofield_r.png', dpi=300)
    plt.close()

    # legend
    pcm = plt.pcolormesh(X, T[I_R_max[0]:Nt - (I_L_max[0] + 3)], Z_sim_i[(I_L_max[0] + 3):Nt - I_R_max[0], :],
                         cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm)
    cbar.set_label('$A\'_I(x,t\')$', rotation=0, size=25, labelpad=-27, y=1.11)

    # put the major ticks at the middle of each cell
    plt.xlim([X[0], X[-1]])
    plt.xlabel('$x$', size='25')
    plt.xticks(fontsize=15)

    plt.ylim([0, tf - ti])
    plt.ylabel('$t\'$', size='25')
    plt.yticks(fontsize=15)

    plt.grid(linestyle='--', alpha=0.2, color='k')
    cbar.ax.tick_params(labelsize=15)

    # labels
    plt.tight_layout()
    plt.savefig(directory + '/Ofield_i.png', dpi=300)
    plt.close()

    # --------------- #

    # legend
    pcm = plt.pcolormesh(X, T[I_R_max[0]:Nt - (I_L_max[0] + 3)], np.abs(Z_r[I_R_max[0]:Nt - (I_L_max[0] +3), :] - Z_sim_r[(I_L_max[0] + 3):Nt - I_R_max[0], :]) / np.mean(np.amax(Z_r)), cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm)
    cbar.set_label('$\delta/A_R^{\\textrm{max}}$', rotation=0, size=25, labelpad=-27, y=1.11)

    # put the major ticks at the middle of each cell
    plt.xlim([X[0], X[-1]])
    plt.xlabel('$x$', size='25')
    plt.xticks(fontsize=15)

    plt.ylim([0, tf - ti])
    plt.ylabel('$t$', size='25')
    plt.yticks(fontsize=15)

    plt.grid(linestyle='--', alpha=0.2, color='k')
    cbar.ax.tick_params(labelsize=15)

    # labels
    plt.tight_layout()
    plt.savefig(directory + '/difference_r.png', dpi=300)
    plt.close()

    # legend
    pcm = plt.pcolormesh(X, T[I_R_max[0]:Nt - (I_L_max[0] + 3)], np.abs(Z_i[I_R_max[0]:Nt - (I_L_max[0] +3), :] - Z_sim_i[(I_L_max[0] + 3):Nt - I_R_max[0], :]) / np.mean(np.amax(Z_i)), cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm)
    cbar.set_label('$\delta/A_I^{\\textrm{max}}$', rotation=0, size=25, labelpad=-27, y=1.11)

    # put the major ticks at the middle of each cell
    plt.xlim([X[0], X[-1]])
    plt.xlabel('$x$', size='25')
    plt.xticks(fontsize=15)

    plt.ylim([0, tf - ti])
    plt.ylabel('$t$', size='25')
    plt.yticks(fontsize=15)

    plt.grid(linestyle='--', alpha=0.2, color='k')
    cbar.ax.tick_params(labelsize=15)

    # labels
    plt.tight_layout()
    plt.savefig(directory + '/difference_i.png', dpi=300)
    plt.close()