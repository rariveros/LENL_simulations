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
    data_out = "D:/mnustes_science/simulation_data/FD/PDNLS_extended_PT/oscillations_ddist_02/data_out"
    index_01 = []
    index_02 = []
    index_03_real = []
    index_03_imag = []
    freqs = []
    phases = []
    parameters = []


    def max_finder_local(Z, Nt, dt):
        D = sparse_D_neumann_4order(Nt, dt)
        DD = sparse_DD_neumann(Nt, dt)
        Z = filtro_array(10, Z)

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

    for directory in directories:
        print("#############   " + directory + "   #############")
        Z_r = np.loadtxt(working_directory + '/' + directory + '/field_real.txt', delimiter=',')
        Z_i = np.loadtxt(working_directory + '/' + directory + '/field_img.txt', delimiter=',')
        X = np.loadtxt(working_directory + '/' + directory + '/X.txt', delimiter=',')
        T = np.loadtxt(working_directory + '/' + directory + '/T.txt', delimiter=',')
        params = np.loadtxt(working_directory + '/' + directory + '/parameters.txt', delimiter=',')
        center_L = - 0.5 * params[3]
        center_R = 0.5 * params[3]

        Nx = len(X)
        Nt = len(T)

        I_L = np.where(X == int(center_L))[0][0]
        I_R = np.where(X == int(center_R))[0][0]

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
        tau = np.arange(-T[-1], T[-1] + dt, dt) # np.arange(-T[-1], T[-1]+ dt, dt)

        ########### SYMMETRY ###########
        middle_tau = int(len(tau) / 2)
        middle_x = int(len(X) / 2)
        S = np.transpose(CCF)
        S_minus = np.transpose(CCF_flipped)
        S_Lx = S[middle_tau, :middle_x-1]
        S_Rx = S[middle_tau, middle_x:-1]
        X_L = X[:middle_x-1]
        X_R = X[middle_x:-1]

        S_L_integred = float(integrate.simpson(S_Lx, X_L))
        S_R_integred = float(integrate.simpson(S_Rx, X_R))
        delta_S = S_Lx - np.flip(S_Rx)

        int_parameter = np.abs((S_L_integred - S_R_integred) / (S_L_integred + S_R_integred))   #va de cero a uno
        desvesta_parameter = np.mean(delta_S ** 2) #no normalizado
        index_01.append(int_parameter)
        index_02.append(desvesta_parameter)
        parameters.append(params[3]) #6:sigma 5:nu 3:dist

        ########### PERIODICITY ###########
        Nt = len(tau)
        dt = tau[1] - tau[0]

        SL_points_max, tau_L_points_max = max_finder_local(S[:, I_L], Nt, dt)
        SR_points_max, tau_R_points_max = max_finder_local(S[:, I_R], Nt, dt)
        SminusL_points_max, tau_min_L_points_max = max_finder_local(S_minus[:, I_L], Nt, dt)
        SminusR_points_max, tau_min_R_points_max = max_finder_local(S_minus[:, I_R], Nt, dt)

        plt.plot(tau, S[:, I_L], color="k", zorder=0)
        plt.scatter(tau_L_points_max, SL_points_max)
        plt.xlabel('$t$', size='20')
        plt.xticks(fontsize=15)
        plt.xlim([tau[0], tau[-1]])
        plt.ylabel('$S(x_L, \\tau)$', size='20')
        plt.yticks(fontsize=15)
        plt.ylim([0, 1.1 * np.amax(S[:, I_L])])
        plt.grid(alpha=0.2)
        plt.tight_layout()
        #plt.show()
        # plt.savefig(analysis_directory + '/S+_L_periods.png', dpi=300)
        plt.close()

        #SL_points_max, tau_L_points_max, I_max_SL = max_finder(S[:, I_L], T, Nt, dt)
        #SR_points_max, tau_R_points_max, I_max_SR = max_finder(S[:, I_R], T, Nt, dt)
        #SminusL_points_max, tau_min_L_points_max, I_max_SminusL = max_finder(S_minus[:, I_L], T, Nt, dt)
        #SminusR_points_max, tau_min_R_points_max, I_max_SminusR = max_finder(S_minus[:, I_R], T, Nt, dt)

        dT_L = np.diff(tau_L_points_max[1:-1])
        dT_R = np.diff(tau_R_points_max[1:-1])

        period_L = np.mean(dT_L)
        period_L_std = np.std(dT_L)
        period_R = np.mean(dT_R)
        period_R_std = np.std(dT_R)
        freq_L = 1 / period_L
        freq_L_std = period_L_std / period_L ** 2
        freq_R = 1 / period_R
        freq_R_std = period_R_std / period_R ** 2

        phase_L = np.mean((tau_min_R_points_max[1:-1] % period_L) / period_L) * 2
        phase_L_std = np.std((tau_min_R_points_max[1:-1] % period_L) / period_L) * 2
        phase_R = np.mean((tau_min_R_points_max[1:-1] % period_R) / period_R) * 2
        phase_R_std = np.std((tau_min_R_points_max[1:-1] % period_L) / period_L) * 2

        freqs.append([freq_L, freq_L_std, freq_R, freq_R_std])
        phases.append([phase_L, phase_L_std, phase_R, phase_R_std])

        ########### SYMMETRY II ###########

        if np.isnan(period_L):
            period_index = 0
            index_03_real.append(0)
            index_03_imag.append(0)
        else:
            period_index = int(period_L / dt)
            print(period_index)
            delta = 0
            Z_delta = Z[:-(int(0.5 * period_index) + delta), :]
            Z_conj_delta = Z_conj[:-(int(0.5 * period_index) + delta), :]

            Z_sim_r = np.flip(Z_r, axis=1)
            Z_sim_i = np.flip(Z_i, axis=1)
            Z_sim_delta = Z_sim_r + 1j * Z_sim_i
            Z_sim_conj_delta = Z_sim_r - 1j * Z_sim_i
            Z_sim_delta = Z_sim_delta[(int(0.5 * period_index) + delta):, :]
            Z_sim_conj_delta = Z_sim_conj_delta[(int(0.5 * period_index) + delta):, :]

            T_delta = T[:Nt - (int(0.5 * period_index) + delta)]

            integral_01 = integrate.simpson(Z_conj_delta * 1j * Z_sim_delta, X)
            integral_02 = integrate.simpson(Z_conj_delta * Z_delta, X)
            sym_02_R = np.mean(np.real(integral_01 / integral_02))
            sym_02_I = np.mean(np.imag(integral_01 / integral_02))
            print([sym_02_R, sym_02_I])
            index_03_real.append(np.real(sym_02_R))
            index_03_imag.append(np.imag(sym_02_I))

    freqs = np.array(freqs)
    phases = np.array(phases)

    ############ SAVING DATA ################
    if not os.path.exists(data_out):
        os.makedirs(data_out)

    np.savetxt(data_out + '/index_01.txt', index_01, delimiter=',')
    np.savetxt(data_out + '/index_02.txt', index_02, delimiter=',')
    np.savetxt(data_out + '/index_03_real.txt', index_03_real, delimiter=',')
    np.savetxt(data_out + '/index_03_imag.txt', index_03_imag, delimiter=',')
    np.savetxt(data_out + '/freqs.txt', freqs, delimiter=',')
    np.savetxt(data_out + '/phases.txt', phases, delimiter=',')
    np.savetxt(data_out + '/parameters.txt', parameters, delimiter=',')