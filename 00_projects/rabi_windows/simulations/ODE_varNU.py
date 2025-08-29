from functions import *
from back_process import *
from time_integrators import *
from numpy.fft import fft, ifft

def gaussian(x, A, x_0, sigma):
    return A * np.exp(- (x - x_0) ** 2 / (2 * sigma ** 2))


def resample(old_arrays, old_x, N_resample):
    new_arrays = []
    for i in range(len(old_arrays)):
        [new_array_i, new_x] = signal.resample(old_arrays[i], int(len(old_x) * (N_resample)), old_x)
        new_arrays.append(new_array_i)
    return new_arrays, new_x


if __name__ == '__main__':

    eq = 'PT_dimer'
    t_rate = 1
    dt = 0.4
    T = 5000
    gamma_str = "0.280"
    sigme_str = "3.000"
    disco = 'D:/'
    initial_dir_data = str(disco) + 'Users/mnustes_science/PT_fluids/mnustes_science/simulation_data'
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    directories = [name for name in os.listdir(working_directory) if
                   os.path.isdir(os.path.join(working_directory, name))]
    for directory in directories:
        Z_r_00 = np.loadtxt(working_directory + '/' + directory + '/sigma=' + sigme_str +'/gamma=' + gamma_str + '/phi00/field_real_0.txt', delimiter=',')
        Z_i_00 = np.loadtxt(working_directory + '/' + directory + '/sigma=' + sigme_str +'/gamma=' + gamma_str + '/phi00//field_img_0.txt', delimiter=',')
        params_00 = np.loadtxt(working_directory + '/' + directory + '/sigma=' + sigme_str +'/gamma=' + gamma_str + '/phi00/parameters.txt', delimiter=',')

        Z_r_01 = np.loadtxt(working_directory + '/' + directory + '/sigma=' + sigme_str +'/gamma=' + gamma_str + '/phi01/field_real_0.txt', delimiter=',')
        Z_i_01 = np.loadtxt(working_directory + '/' + directory + '/sigma=' + sigme_str +'/gamma=' + gamma_str + '/phi01/field_img_0.txt', delimiter=',')
        x_grid = np.loadtxt(working_directory + '/' + directory + '/sigma=' + sigme_str +'/gamma=' + gamma_str + '/phi01/X.txt', delimiter=',')
        params_01 = np.loadtxt(working_directory + '/' + directory + '/sigma=' + sigme_str +'/gamma=' + gamma_str + '/phi01/parameters.txt', delimiter=',')

        dx = x_grid[1] - x_grid[0]
        #[Z_r_00, Z_i_00, Z_r_01, Z_i_01] = [Z_r_00[-1, :], Z_i_00[-1, :], Z_r_01[-1, :], Z_i_01[-1, :]]

        distances = np.arange(1, 40, 0.25)
        PHI_Rs = []
        PHI_Ls = []
        times = []
        UR_Rs = []
        UR_Is = []
        UL_Rs = []
        UL_Is = []
        U1_Rs = []
        U1_Is = []
        U2_Rs = []
        U2_Is = []

        SIGMAS = []
        PIS = []
        DELTAS = []
        DELTAS1 = []
        DELTAS2 = []
        nu = params_01[4]
        for d in distances:
            print("############ nu = " + str(nu) + " ############")
            print("### d = " + str(d) + " ###")
            X_L = - d / 2
            X_R = + d / 2
            J_L = np.argmin(np.abs(x_grid - X_L))
            J_R = np.argmin(np.abs(x_grid - X_R))
            J_center = np.argmin(np.abs(x_grid))
            Delta_J_L = J_center - J_L
            Delta_J_R = J_R - J_center
            PHI_L = np.append(Z_r_00[Delta_J_L:], np.zeros(Delta_J_L)) + 1j * np.append(Z_i_00[Delta_J_L:], np.zeros(Delta_J_L))
            PHI_R = np.append(np.zeros(Delta_J_R), Z_r_01[:-Delta_J_R]) + 1j * np.append(np.zeros(Delta_J_R), Z_i_01[:-Delta_J_R])
            PHI_L = 1j * PHI_L
            PHI_R = 1j * PHI_R

            PHI = PHI_L + PHI_R
            phi_01 = PHI_L - 1j * PHI_R
            phi_02 = PHI_L + 1j * PHI_R

            [alpha, beta, gamma_0, mu, nu, sigma, phi] = params_00
            gamma = -gamma_0 * (np.exp(- (x_grid - d / 2) ** 2 / (2 * sigma ** 2)) - np.exp(- (x_grid + d / 2) ** 2 / (2 * sigma ** 2)))

            Gamma_0 = integrate.simpson(phi_01 * phi_02, x_grid)
            Gamma_1 = integrate.simpson(phi_01 ** 2, x_grid)
            Gamma_2 = integrate.simpson(phi_02 ** 2, x_grid)

            ddx_phi_01 = alpha * np.append(0, np.append(np.diff(np.diff(phi_01)) / (dx ** 2), 0))
            ddx_phi_02 = alpha * np.append(0, np.append(np.diff(np.diff(phi_02)) / (dx ** 2), 0))

            Sigma_11 = integrate.simpson(phi_01 * ddx_phi_01, x_grid) / Gamma_1
            Sigma_12 = integrate.simpson(phi_01 * ddx_phi_02, x_grid) / Gamma_1
            Sigma_21 = integrate.simpson(phi_02 * ddx_phi_01, x_grid) / Gamma_2
            Sigma_22 = integrate.simpson(phi_02 * ddx_phi_02, x_grid) / Gamma_2

            Pi_11 = integrate.simpson(phi_01 * np.conjugate(phi_01) * gamma, x_grid) / Gamma_1
            Pi_12 = integrate.simpson(phi_01 * np.conjugate(phi_02) * gamma, x_grid) / Gamma_1
            Pi_21 = integrate.simpson(phi_02 * np.conjugate(phi_01) * gamma, x_grid) / Gamma_2
            Pi_22 = integrate.simpson(phi_02 * np.conjugate(phi_02) * gamma, x_grid) / Gamma_2

            Delta_11 = (1) * integrate.simpson(phi_01 * np.abs(phi_02) ** 2 * phi_02, x_grid) / Gamma_1
            Delta_21 = (2) * integrate.simpson(phi_01 * np.abs(phi_02) ** 2 * phi_01, x_grid) / Gamma_1
            Delta_31 = (2) * integrate.simpson(phi_01 * np.abs(phi_01) ** 2 * phi_02, x_grid) / Gamma_1
            Delta_41 = (1) * integrate.simpson(phi_01 * np.abs(phi_01) ** 2 * phi_01, x_grid) / Gamma_1
            Delta_51 = (1) * integrate.simpson(phi_01 * phi_02 ** 2 * np.conjugate(phi_01), x_grid) / Gamma_1
            Delta_61 = (1) * integrate.simpson(phi_01 * phi_01 ** 2 * np.conjugate(phi_02), x_grid) / Gamma_1

            Delta_12 = (1) * integrate.simpson(phi_02 * np.abs(phi_02) ** 2 * phi_02, x_grid) / Gamma_2
            Delta_22 = (2) * integrate.simpson(phi_02 * np.abs(phi_02) ** 2 * phi_01, x_grid) / Gamma_2
            Delta_32 = (2) * integrate.simpson(phi_02 * np.abs(phi_01) ** 2 * phi_02, x_grid) / Gamma_2
            Delta_42 = (1) * integrate.simpson(phi_02 * np.abs(phi_01) ** 2 * phi_01, x_grid) / Gamma_2
            Delta_52 = (1) * integrate.simpson(phi_02 * phi_02 ** 2 * np.conjugate(phi_01), x_grid) / Gamma_2
            Delta_62 = (1) * integrate.simpson(phi_02 * phi_01 ** 2 * np.conjugate(phi_02), x_grid) / Gamma_2

            dJ = np.array([[- mu + np.imag(Sigma_11) + np.real(Pi_11), - nu - np.real(Sigma_11) + np.imag(Pi_11)  , np.imag(Sigma_21) + np.real(Sigma_21), -np.real(Sigma_21) + np.imag(Pi_21)],
                           [+ nu + np.real(Sigma_11) + np.imag(Pi_11), - mu + np.imag(Sigma_11) - np.real(Pi_11), np.real(Sigma_21) + np.imag(Pi_21), np.real(Sigma_21) - np.real(Pi_21)],
                           [np.imag(Sigma_12) + np.real(Pi_12)       ,  -np.real(Sigma_12) + np.imag(Pi_12)      , - mu + np.imag(Sigma_22) + np.real(Pi_22), -nu - np.real(Sigma_22) + np.imag(Pi_22)],
                           [np.real(Sigma_12) + np.imag(Pi_12), np.imag(Sigma_12) - np.real(Pi_12)             , +nu + np.real(Sigma_22) + np.real(Pi_22), - mu + np.imag(Sigma_22) + np.imag(Pi_22)]])
            eigenvalues, eigenvectors = np.linalg.eig(dJ)
            eig_R = np.real(eigenvalues)
            eig_I = np.imag(eigenvalues)
            for i in [0, 1, 2, 3]:
                if i == 0 and d == distances[0]:
                    plt.scatter(d, eig_R[i], c="b", s=10, label="$\lambda_R$")
                    plt.scatter(d, eig_I[i], c="r", s=10, label="$\lambda_I$")
                else:
                    plt.scatter(d, eig_R[i], c="b", s=10)
                    plt.scatter(d, eig_I[i], c="r", s=10)

            sigma = [Sigma_11, Sigma_12, Sigma_21, Sigma_22]
            pi = [Pi_11, Pi_12, Pi_21, Pi_22]
            delta1 = [Delta_11, Delta_21, Delta_31, Delta_41, Delta_51, Delta_61]
            delta2 = [Delta_12, Delta_22, Delta_32, Delta_42, Delta_52, Delta_62]
            delta = [Delta_11, Delta_21, Delta_31, Delta_41, Delta_51, Delta_61, Delta_12, Delta_22, Delta_32, Delta_42, Delta_52, Delta_62]

            SIGMAS.append(sigma)
            PIS.append(pi)
            DELTAS1.append(delta1)
            DELTAS2.append(delta2)

            parameters = [alpha, beta, mu, nu, sigma, pi, delta]

            # Definición de la grilla
            [tmin, tmax, dt] = [0, T, dt]
            t_grid = np.arange(tmin, tmax + dt, dt)
            T = tmax
            Nt = t_grid.shape[0]

            # Initial Conditions Pattern
            U_1_init = 0.1 * (np.random.rand(1) + 1j * np.random.rand(1))
            U_2_init = 0.1 * (np.random.rand(1) + 1j * np.random.rand(1))
            operators = [0]

            # Empaquetamiento de parametros, campos y derivadas para integración
            fields_init = [U_1_init, U_2_init]
            grids = [t_grid, 0, 0]

            # Midiendo tiempo inicial
            now = datetime.datetime.now()
            print('Hora de Inicio: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
            time_init = time.time()

            final_fields, fields_history, time_grid = RK4_FD(eq, fields_init, parameters, grids, dt, Nt, operators, t_rate) #INTEGRACION EN EL TIEMPO

            now = datetime.datetime.now()
            print('Hora de Término: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
            time_fin = time.time()
            print(str(time_fin - time_init) + ' seg')

            # Reobteniendo campos
            U1_light = np.array(fields_history)[:, 0]
            U2_light = np.array(fields_history)[:, 1]
            t_light = time_grid

            UL = (U1_light + U2_light)
            UR = - 1j * (U1_light - U2_light)
            U1 = U1_light
            U2 = U2_light

            U1_flatten = []
            U2_flatten = []
            UL_flatten = []
            UR_flatten = []
            for i in range(len(UL)):
                UL_flatten.append(UL[i, 0])
                UR_flatten.append(UR[i, 0])
                U1_flatten.append(U1[i, 0])
                U2_flatten.append(U2[i, 0])
            UL = UL_flatten
            UR = UR_flatten
            U1 = U1_flatten
            U2 = U2_flatten
            PHI_Rs.append(PHI_R)
            PHI_Ls.append(PHI_L)
            UR_Rs.append(np.real(UR))
            UR_Is.append(np.imag(UR))
            UL_Rs.append(np.real(UL))
            UL_Is.append(np.imag(UL))
            U1_Rs.append(np.real(U1))
            U1_Is.append(np.imag(U1))
            U2_Rs.append(np.real(U2))
            U2_Is.append(np.imag(U2))

        save_directory = working_directory + '/' + directory + '/sigma=' + sigme_str +'/gamma=' + gamma_str + "/analysis"
        if not os.path.exists(save_directory):
            os.makedirs(save_directory)

        plt.grid(alpha=0.2)
        plt.xticks(fontsize=15)
        plt.xlabel("$d$", fontsize=20)
        plt.yticks(fontsize=15)
        plt.ylabel("$\lambda_R, \lambda_I$", fontsize=20)
        plt.legend(fontsize=15)
        plt.tight_layout()
        plt.savefig(save_directory + "/eigen.png", dpi=200)
        plt.close()

        labels_01 = ["$\Sigma_{++}$", "$\Sigma_{+-}$", "$\Sigma_{-+}$", "$\Sigma_{--}$"]
        labels_02 = ["$\Pi_{++}$", "$\Pi_{+-}$", "$\Pi_{-+}$", "$\Pi_{--}$"]
        #[Delta_11, Delta_21, Delta_31, Delta_41, Delta_51, Delta_61, Delta_12, Delta_22, Delta_32, Delta_42, Delta_52, Delta_62]
        labels_03 = ["$\Delta_1^{+}$", "$\Delta_2^{+}$", "$\Delta_3^{+}$", "$\Delta_4^{+}$", "$\Delta_5^{+}$", "$\Delta_6^{+}$"]
        labels_04 = ["$\Delta_1^{-}$", "$\Delta_2^{-}$", "$\Delta_3^{-}$", "$\Delta_4^{-}$", "$\Delta_5^{-}$", "$\Delta_6^{-}$"]

        fig, ((ax11, ax12), (ax21, ax22), (ax31, ax32), (ax41, ax42)) = plt.subplots(4, 2, figsize=(7, 12))
        for i in range(4):
            ax11.plot(distances, np.real(np.array(SIGMAS)[:, i]), label=labels_01[i])
            #ax11.legend(fontsize=10, loc="upper right")
            ax12.plot(distances, np.imag(np.array(SIGMAS)[:, i]), label=labels_01[i])
            ax12.legend(fontsize=10, loc="upper right")

        for i in range(4):
            ax21.plot(distances, np.real(np.array(PIS)[:, i]), label=labels_02[i])
            #ax21.plot(distances, np.imag(np.array(PIS)[:, i]), label=labels_02[i], ls="--")
            #ax21.legend(fontsize=10, loc="upper right")
            ax22.plot(distances, np.imag(np.array(PIS)[:, i]), label=labels_02[i])
            ax22.legend(fontsize=10, loc="upper right")

        for i in range(6):
            ax31.plot(distances, np.real(np.array(DELTAS1)[:, i]), label=labels_03[i])
            #ax31.legend(fontsize=10, loc="upper right")
            ax32.plot(distances, np.imag(np.array(DELTAS1)[:, i]), label=labels_03[i])
            ax32.legend(fontsize=10, loc="upper right")

        for i in range(6):
            ax41.plot(distances, np.real(np.array(DELTAS2)[:, i]), label=labels_04[i])
            #ax41.legend(fontsize=10, loc="upper right")
            ax42.plot(distances, np.imag(np.array(DELTAS2)[:, i]), label=labels_04[i])
            ax42.legend(fontsize=10, loc="upper right")
        plt.savefig("coefficients.png", dpi=300)
        plt.close()


        np.savetxt(save_directory + '/ansatz_right.txt', PHI_Rs, delimiter=',')
        np.savetxt(save_directory + '/ansatz_left.txt', PHI_Ls, delimiter=',')
        np.savetxt(save_directory + '/t_grid.txt', t_light, delimiter=',')
        np.savetxt(save_directory + '/x_grid.txt', x_grid, delimiter=',')
        np.savetxt(save_directory + '/UR_Rs.txt', np.array(UR_Rs), delimiter=',')
        np.savetxt(save_directory + '/UR_Is.txt',np.array(UR_Is), delimiter=',')
        np.savetxt(save_directory + '/UL_Rs.txt', np.array(UL_Rs), delimiter=',')
        np.savetxt(save_directory + '/UL_Is.txt', np.array(UL_Is), delimiter=',')
        np.savetxt(save_directory + '/U1_Rs.txt', np.array(U1_Rs), delimiter=',')
        np.savetxt(save_directory + '/U1_Is.txt',np.array(U1_Is), delimiter=',')
        np.savetxt(save_directory + '/U2_Rs.txt', np.array(U2_Rs), delimiter=',')
        np.savetxt(save_directory + '/U2_Is.txt', np.array(U2_Is), delimiter=',')
        np.savetxt(save_directory + '/dists.txt', distances, delimiter=',')
