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
    data_out = "D:/mnustes_science/simulation_data/FD/PDNLS_extended_PT/dgamma/data_out"
    PL = []
    AL = []
    PR = []
    AR = []
    AL_error = []
    AR_error = []
    parameters = []

    def max_finder(Z, T):
        dt = T[1] - T[0]
        Nt = len(T)
        D = sparse_D_neumann_4order(Nt, dt)
        DD = sparse_DD_neumann(Nt, dt)
        D1_Z = D.dot(np.transpose(Z))
        D2_Z = DD.dot(np.transpose(Z))

        tau_L_points_max = []
        Z_points_max = []
        for i in range(Nt):
            if np.sign(D1_Z[i]) != np.sign(D1_Z[i - 1]) and D2_Z[i] < 0 and i != 0 and np.sign(D1_Z[i]) + np.sign(
                    D1_Z[i - 1]) != 1:
                tau_L_points_max.append(T[i])
                Z_points_max.append(Z[i])
        Z_points_max = np.array(Z_points_max)
        tau_L_points_max = np.array(tau_L_points_max)
        return Z_points_max, tau_L_points_max

    def min_finder(Z, T):
        dt = T[1] - T[0]
        Nt = len(T)
        D = sparse_D_neumann_4order(Nt, dt)
        DD = sparse_DD_neumann(Nt, dt)
        D1_Z = D.dot(np.transpose(Z))
        D2_Z = DD.dot(np.transpose(Z))

        tau_L_points_min = []
        Z_points_min = []
        for i in range(Nt):
            if np.sign(D1_Z[i]) != np.sign(D1_Z[i - 1]) and D2_Z[i] > 0 and i != 0 and np.sign(D1_Z[i]) + np.sign(
                    D1_Z[i - 1]) != -1:
                tau_L_points_min.append(T[i])
                Z_points_min.append(Z[i])
        Z_points_min = np.array(Z_points_min)
        tau_L_points_min = np.array(tau_L_points_min)
        return Z_points_min, tau_L_points_min

    for directory in directories:
        print("#############   " + directory + "   #############")
        Z_r = np.loadtxt(working_directory + '/' + directory + '/dist=20.000/field_real.txt', delimiter=',')
        Z_i = np.loadtxt(working_directory + '/' + directory + '/dist=20.000/field_img.txt', delimiter=',')
        X = np.loadtxt(working_directory + '/' + directory + '/dist=20.000/X.txt', delimiter=',')
        T = np.loadtxt(working_directory + '/' + directory + '/dist=20.000/T.txt', delimiter=',')
        params = np.loadtxt(working_directory + '/' + directory + '/dist=20.000/parameters.txt', delimiter=',')
        center_L = - 0.5 * params[3]
        center_R = 0.5 * params[3]

        Nx = len(X)
        Nt = len(T)

        I_L = np.where(X == center_L)[0][0]
        I_R = np.where(X == center_R)[0][0]

        ti, tf = -400,  -1 #int(0.7 * Nt), Nt
        Z_r = Z_r[ti:tf, :]
        Z_i = Z_i[ti:tf, :]
        T = T[ti:tf] - T[ti]
        dt = T[1] - T[0]

        Nx = len(X)
        Nt = len(T)

        Z = Z_r + 1j * Z_i
        Z_conj = Z_r - 1j * Z_i
        Z_modulo = np.absolute(Z)

        Z_L_points_min, T_L_points_min = min_finder(Z_modulo[:, I_L], T)
        Z_L_points_max, T_L_points_max = max_finder(Z_modulo[:, I_L], T)
        Z_R_points_min, T_R_points_min = min_finder(Z_modulo[:, I_R], T)
        Z_R_points_max, T_R_points_max = max_finder(Z_modulo[:, I_R], T)

        #plt.plot(T, Z_modulo[:, I_L])
        #plt.scatter(T_L_points_min, Z_L_points_min)
        #plt.title("$\gamma=" + str(params[2]) + "\  \ N_t=" + str(Nt) + "$")
        #plt.show()

        AL_i = 0.5 * (np.mean(Z_L_points_max) - np.mean(Z_L_points_min))
        AL_error_i = 0.5 * np.sqrt(np.std(Z_L_points_max) ** 2 * np.mean(Z_L_points_max) + np.std(Z_L_points_min) ** 2 * np.mean(Z_L_points_min))
        AR_i = 0.5 * (np.mean(Z_R_points_max) - np.mean(Z_R_points_min))
        AR_error_i = 0.5 * np.sqrt(
            np.std(Z_R_points_max) ** 2 * np.mean(Z_R_points_max) + np.std(Z_R_points_min) ** 2 * np.mean(
                Z_R_points_min))

        AL.append([AL_i, AL_error_i])
        AR.append([AR_i, AR_error_i])
        parameters.append(params[2])
    AL = np.array(AL)
    AR = np.array(AR)
    if not os.path.exists(data_out):
        os.makedirs(data_out)

    fig, ax = plt.subplots()
    ax.errorbar(parameters, AL[:, 0], AL[:, 1], marker='o', ls='', ecolor="k", mec='black',
                color=(216 / 255, 135 / 255, 25 / 255), label="$x_L$")
    plt.xlabel('$\sigma_i$', size='25')
    plt.ylabel('$\phi$', size='25', rotation=0)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    #plt.xlim(2.9, 6.5)
    #plt.ylim(0, 2)
    #ax.yaxis.set_major_formatter(FormatStrFormatter('$%g \pi$'))
    #ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base=0.25))
    plt.grid(linestyle='--', alpha=0.5)
    #plt.legend(loc="upper right", fontsize=20)
    plt.tight_layout()
    plt.show()
    #plt.savefig(directory + '//phases.png', dpi=300)
    plt.close()

    np.savetxt(data_out + '/AR_up.txt', AR, delimiter=',')
    np.savetxt(data_out + '/AL_up.txt', AL, delimiter=',')
    np.savetxt(data_out + '/params_up.txt', parameters, delimiter=',')

