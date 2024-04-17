from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    disco = 'D:/'
    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FDa'
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    directories = [name for name in os.listdir(working_directory) if os.path.isdir(os.path.join(working_directory, name))]
    data_out = "C:/mnustes_science/simulation_data/FD/PDNLS_extended_PT/extras/dimensional/analysis"
    index_01 = []
    index_02 = []
    parameters = []

    for directory in directories:
        print("#############   " + directory + "   #############")
        Z_r = np.loadtxt(working_directory + '/' + directory + 'field_real.txt', delimiter=',')
        Z_i = np.loadtxt(working_directory + '/' + directory + '/field_img.txt', delimiter=',')
        X = np.loadtxt(working_directory + '/' + directory + '/X.txt', delimiter=',')
        T = np.loadtxt(working_directory + '/' + directory + '/T.txt', delimiter=',')
        params = np.loadtxt(working_directory + '/' + directory + '/parameters.txt', delimiter=',')
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

        delta = 0
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

        index_01.append(np.mean(np.real(integral_01 / integral_02)))
        index_02.append(np.mean(np.imag(integral_01 / integral_02)))
        parameters.append(params[3])

    if not os.path.exists(data_out):
        os.makedirs(data_out)

    plt.scatter(parameters, index_01)
    plt.scatter(parameters, index_02)
    plt.xlabel('$\\nu$', size='20')
    plt.ylabel('$\\xi^{1}_s$', size='20')
    #plt.ylim(0.05, 0.125)
    plt.grid(linestyle='--', alpha=0.5)
    plt.show()
    #plt.savefig(data_out + '//index_01.png', dpi=300)
    plt.close()

