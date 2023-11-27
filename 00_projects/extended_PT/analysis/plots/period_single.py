from functions import *
from back_process import *
from time_integrators import *


if __name__ == '__main__':
    disco = 'D:/'
    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FDa'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    def max_finder_local(Z, Nt, dt, tau):
        D = sparse_D_neumann_4order(Nt, dt)
        DD = sparse_DD_neumann(Nt, dt)
        #Z = filtro_array(10, Z)

        D1_Z = D.dot(np.transpose(Z))
        D2_Z = DD.dot(np.transpose(Z))

        tau_L_points_max = []
        Z_points_max = []
        for i in range(Nt):
            if np.sign(D1_Z[i]) != np.sign(D1_Z[i - 1]) and D2_Z[i] < 0 and i != 0 and np.sign(D1_Z[i]) + np.sign(
                    D1_Z[i - 1]) != 1:
                tau_L_points_max.append(tau[i])
                Z_points_max.append(Z[i])
        Z_points_max = np.array(Z_points_max)
        tau_L_points_max = np.array(tau_L_points_max)
        return Z_points_max, tau_L_points_max

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

    ti, tf = int(0.5 * Nt), Nt
    Z_r = Z_r[ti:tf, :]
    Z_i = Z_i[ti:tf, :]
    T = T[ti:tf] - T[ti]
    dt = T[1] - T[0]

    Nx = len(X)
    Nt = len(T)

    Z = Z_r + 1j * Z_i
    Z_conj = Z_r - 1j * Z_i
    Z_modulo = np.absolute(Z)

    CCF = []
    CCF_flipped = []

    for j in range(Nx):
        CCF_j = np.correlate(Z[:, j], Z_conj[:, j], "full")
        CCF_flipped_j = np.correlate(np.flip(Z, axis=1)[:, j], Z_conj[:, j], "full")
        CCF.append(CCF_j)
        CCF_flipped.append(CCF_flipped_j)
    CCF = np.array(np.abs(CCF))
    CCF_flipped = np.array(np.abs(CCF_flipped))
    tau = np.arange(-T[-1], T[-1] + dt, dt)  # np.arange(-T[-1], T[-1]+ dt, dt)

    ########### SYMMETRY ###########
    middle_tau = int(len(tau) / 2)
    middle_x = int(len(X) / 2)
    S = np.transpose(CCF)

    ########### PERIODICITY ###########
    Nt = len(tau)
    dt = tau[1] - tau[0]
    print(np.shape(tau))
    print(np.shape(S[:, I_L]))

    SL_points_max, tau_L_points_max = max_finder_local(S[:, I_L], Nt-1, dt, tau[:-1])
    SR_points_max, tau_R_points_max = max_finder_local(S[:, I_R], Nt-1, dt, tau[:-1])

    plt.plot(tau[:-1], S[:, I_L], color="k", zorder=1)
    plt.scatter(tau_L_points_max, SL_points_max, facecolors='r', edgecolors='k', s=50, zorder=2)
    plt.vlines(tau_L_points_max, 0, 1.2 * np.amax(S[:, I_L]), color="r", zorder=0, linestyles="--", alpha=0.4)
    plt.xlabel('$\\tau$', size='20')
    plt.xticks(fontsize=15)
    plt.xlim([tau[0], tau[-1]])
    plt.ylabel('$S(x_L, \\tau)$', size='20')
    plt.yticks(fontsize=15)
    plt.ylim([0, 1.1 * np.amax(S[:, I_L])])
    plt.grid(alpha=0.2)
    plt.tight_layout()
    #plt.savefig(directory + '/SL_periods.png', dpi=200)
    plt.close()

    plt.plot(tau[:-1], S[:, I_R], color="k", zorder=0)
    plt.scatter(tau_R_points_max, SR_points_max, facecolors='r', edgecolors='k', s=50)
    plt.vlines(tau_R_points_max, 0, 1.2 * np.amax(S[:, I_R]), color="r", zorder=0, linestyles="--", alpha=0.4)
    plt.xlabel('$\\tau$', size='20')
    plt.xticks(fontsize=15)
    plt.xlim([tau[0], tau[-1]])
    plt.ylabel('$S(x_L, \\tau)$', size='20')
    plt.yticks(fontsize=15)
    plt.ylim([0, 1.1 * np.amax(S[:, I_L])])
    plt.grid(alpha=0.2)
    plt.tight_layout()
    #plt.savefig(directory + '/SR_periods.png', dpi=200)
    plt.close()