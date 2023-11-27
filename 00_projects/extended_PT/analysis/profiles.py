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

    J_L = np.where(X == center_L)[0][0]
    J_R = np.where(X == center_R)[0][0]

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

    Z_mod_L_points_max, tau_R_points_max, I_R = max_finder(Z_modulo[:, J_L], T, Nt, dt)
    Z_mod_R_points_max, tau_L_points_max, I_L = max_finder(Z_modulo[:, J_R], T, Nt, dt)
    if I_R[0] > I_L[0]:
        I_R_max = I_L
        I_L_max = I_R
    elif I_R[0] < I_L[0]:
        I_R_max = I_R
        I_L_max = I_L
    period = int((np.mean(np.diff(I_L_max[1:-1])) + np.mean(np.diff(I_R_max[1:-1])))/2)

    matplotlib.rc('xtick', labelsize=18)
    matplotlib.rc('ytick', labelsize=18)

    fig, (ax1, ax2) = plt.subplots(2)
    ax1.plot(T, Z_modulo[:, J_R], label="$|A(x_R, t)|$", c="b", linewidth=2)
    ax1.plot(T, Z_modulo[:, J_L], label="$|A(x_L, t)|$", c="r", linewidth=2)
    ax1.set_xlabel("$t$", fontsize=25)
    ax1.set_xlim(0, 200)
    ax1.set_ylim(0, 0.11)#1.1 * np.amax(Z_modulo[I_R[0] + 10, :]))
    #ax1.set_ylabel("$Re(A), Im(A)$", fontsize=25)
    ax1.legend(loc="upper right", fontsize=15)
    ax1.grid(alpha=0.3)

    ax2.plot(X, Z_modulo[I_R[0], :], label="$|A(t)|$", c="indigo", linewidth=2, zorder=2)
    ax2.plot(X, Z_modulo[I_R[0] + int(period / 2), :], label="$|A(t + \\tau_{\\textrm{so}}/2)|$", c="darkorange", linewidth=2, zorder=2)
    ax2.vlines([X[J_L], X[J_R]], 0, 1.1 * np.amax(Z_modulo[I_R[0], :]), linewidths=3, colors="r", linestyles="--", alpha=0.6, zorder=0)
    plt.xlabel("$x\ (\\textrm{mm})$", fontsize=25)
    plt.ylabel("$|A|$", fontsize=25)
    ax2.set_xlim([-75, 75])
    ax2.set_ylim(0, 1.1 * np.amax(Z_modulo[I_R[0], :]))
    ax2.legend(loc="upper right", fontsize=15)
    ax2.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig('typeIB_profile.png', dpi=300)
    plt.close()
