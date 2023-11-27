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

    Z_mod_L_points_max, tau_R_points_max, I_R = max_finder(Z_modulo[:, J_L], T, Nt, dt)
    Z_mod_R_points_max, tau_L_points_max, I_L = max_finder(Z_modulo[:, J_R], T, Nt, dt)
    if I_R[0] > I_L[0]:
        I_R_max = I_L
        I_L_max = I_R
    elif I_R[0] < I_L[0]:
        I_R_max = I_R
        I_L_max = I_L
    Z_sim_r = np.flip(Z_r, axis=1)
    Z_sim_i = np.flip(Z_i, axis=1)

    deltas = np.arange(0, 6)
    for delta in deltas:
        Z_delta = Z[I_R_max[0]:Nt - (I_L_max[0] + delta), :]
        Z_conj_delta = Z_conj[I_R_max[0]:Nt - (I_L_max[0] + delta), :]

        Z_sim_delta = Z_sim_r + 1j * Z_sim_i
        Z_sim_conj_delta = Z_sim_r - 1j * Z_sim_i
        Z_sim_delta = Z_sim_delta[(I_L_max[0] + delta):Nt - I_R_max[0], :]
        Z_sim_conj_delta = Z_sim_conj_delta[(I_L_max[0] + delta):Nt - I_R_max[0], :]

        T_delta = T[I_R_max[0]:Nt - (I_L_max[0] + delta)]

        integral_01 = integrate.simpson(Z_conj_delta * 1j * Z_sim_delta, X)
        integral_02 = integrate.simpson(Z_conj_delta * Z_delta, X)

        plt.plot(T_delta, np.real(integral_01 / integral_02), label="$\\tau_i ="+str(I_L_max[0] + delta)+"$")#$\\textrm{Re}\left\{ \zeta \\right\}$")
        #plt.plot(T_delta, np.imag(integral_01 / integral_02), label="$\\textrm{Im}\left\{ \zeta \\right\}$")
    plt.xlim([0, 100])
    plt.xlabel("$t$", fontsize=20)
    plt.ylabel("$\\textrm{Re}(\zeta)$", fontsize=20)
    plt.legend(loc="upper right")
    plt.tight_layout()
    #plt.show()
    plt.close()

    delta = 3
    Z_delta = Z[I_R_max[0]:Nt - (I_L_max[0] + delta), :]
    Z_conj_delta = Z_conj[I_R_max[0]:Nt - (I_L_max[0] + delta), :]

    Z_sim_delta = Z_sim_r + 1j * Z_sim_i
    Z_sim_conj_delta = Z_sim_r - 1j * Z_sim_i
    Z_sim_delta = Z_sim_delta[(I_L_max[0] + delta):Nt - I_R_max[0], :]
    Z_sim_conj_delta = Z_sim_conj_delta[(I_L_max[0] + delta):Nt - I_R_max[0], :]

    T_delta = T[I_R_max[0]:Nt - (I_L_max[0] + delta)]

    integral_01 = integrate.simpson(Z_conj_delta * 1j * Z_sim_delta, X)
    integral_02 = integrate.simpson(Z_conj_delta * Z_delta, X)
    print(np.mean(np.real(integral_01 / integral_02)))
    print(np.std(np.real(integral_01 / integral_02)))
    print(np.mean(np.imag(integral_01 / integral_02)))
    print(np.std(np.imag(integral_01 / integral_02)))

    fig, ax = plt.subplots()
    plt.plot(T_delta, np.real(integral_01 / integral_02), label="$\\textrm{Re}\left( \zeta \\right)$", c="b", linewidth=3)
    plt.plot(T_delta, np.imag(integral_01 / integral_02), label="$\\textrm{Im}\left( \zeta \\right)$", c="r", linewidth=3)
    plt.xlim([T_delta[0], T_delta[-1]])
    plt.xlabel("$t$", fontsize=30)
    plt.ylabel("$\zeta(t)$", fontsize=30)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.legend(loc="center right", fontsize=18)
    plt.grid(alpha=0.3)
    ax.set_aspect(300)
    plt.tight_layout()
    plt.savefig(directory + '/zeta_t.png', dpi=300)
    plt.close()
