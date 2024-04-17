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
        Z_r = np.loadtxt(working_directory + '/' + directory + '/dist=20.000/field_real.txt', delimiter=',')
        Z_i = np.loadtxt(working_directory + '/' + directory + '/dist=20.000/field_img.txt', delimiter=',')
        X = np.loadtxt(working_directory + '/' + directory + '/dist=20.000/X.txt', delimiter=',')
        T = np.loadtxt(working_directory + '/' + directory + '/dist=20.000/T.txt', delimiter=',')
        params = np.loadtxt(working_directory + '/' + directory + '/dist=20.000/parameters.txt', delimiter=',')

        dt = T[1] - T[0]
        Z = Z_r + 1j * Z_i
        Z_modulo = np.absolute(Z)

        Nx = len(X)
        Nt = len(T)

        ti, tf = int(0.2 * Nt), Nt
        Z_r = Z_r[ti:tf, :]
        Z_i = Z_i[ti:tf, :]
        Z_modulo = Z_modulo[ti:tf, :]
        T = T[ti:tf] - T[ti]
        dt = T[1] - T[0]

        Nx = len(X)
        Nt = len(T)
        Z_complex = Z_r + 1j * Z_i

        ###

        Z_mod_L = Z_modulo[:, 0:int(Nx / 2)]
        Z_mod_R = Z_modulo[:, int(Nx / 2) + 1:]
        X_L = X[0:int(Nx / 2)]
        X_R = X[int(Nx / 2) + 1:]

        N_L = integrate.simpson(Z_mod_L ** 2, X_L)
        N_R = integrate.simpson(Z_mod_R ** 2, X_R)

        CCF = np.correlate(N_L - np.mean(N_L), N_R - np.mean(N_R), "full")
        tau = np.arange(-T[-1], T[-1], dt)
        Ntau = len(tau)
        dtau = tau[1] - tau[0]

        CCF_max, tau_max, CCF_I = max_finder(CCF[:-1], tau, Ntau, dtau)
        #plt.plot(tau, CCF[:-1])
        #plt.show()
        #plt.close()

        tau_R = []
        maxval = np.amax(CCF)
        for i in range(len(tau_max)):
            if CCF_max[i] > 0.5 * maxval:
                tau_R.append(tau_max[i])
        #if len(tau_R) == 2:
        #    period = 0
        #    period_std = 0
        #else:
        period = np.mean(np.diff(tau_R))
        period_std = np.std(np.diff(tau_R))
        print(period)
        t_dist = np.abs(period - T)
        I_period = np.argmin(t_dist)

        Z_conjugated = np.conjugate(Z_complex[:-int(I_period/2), :])
        Z_transformed = 1j * np.flip(Z_complex[int(I_period/2):, :], axis=1)

        integral_01 = integrate.simpson(Z_conjugated * Z_transformed, X)
        integral_02 = integrate.simpson(np.abs(Z_conjugated) ** 2, X)
        print(np.mean(np.real(integral_01 / integral_02)))
        print(np.mean(np.imag(integral_01 / integral_02)))

        index_01.append(np.mean(np.real(integral_01 / integral_02)))
        index_02.append(np.mean(np.imag(integral_01 / integral_02)))
        parameters.append(params[2])

    if not os.path.exists(data_out):
        os.makedirs(data_out)
    np.savetxt(data_out + '/coef_real_gamma.txt', index_01, delimiter=',')
    np.savetxt(data_out + '/coef_imag_gamma.txt', index_02, delimiter=',')
    np.savetxt(data_out + '/coef_gamma.txt', parameters, delimiter=',')

    plt.scatter(parameters, index_01)
    plt.scatter(parameters, index_02)
    plt.xlabel('$\gamma_0$', size='20')
    plt.ylabel('$\\xi_S$', size='20')
    plt.grid(linestyle='--', alpha=0.5)
    plt.savefig(data_out + '/coef_gamma.png', dpi=300)
    plt.close()
