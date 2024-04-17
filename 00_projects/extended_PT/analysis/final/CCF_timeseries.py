import matplotlib.pyplot as plt

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
    DATA = []
    profiles_R = []
    profiles_L = []
    X_domain = []
    m = 0
    for directory in directories:
        m = m + 1
        #directory = "gamma=0.2700"
        print("#############   " + directory + "   #############")
        Z_r = np.loadtxt(working_directory + "/" + directory + '/field_real.txt', delimiter=',')
        Z_i = np.loadtxt(working_directory + "/" + directory + '/field_img.txt', delimiter=',')
        X = np.loadtxt(working_directory + "/" + directory + '/X.txt', delimiter=',')
        T = np.loadtxt(working_directory + "/" + directory + '/T.txt', delimiter=',')
        params = np.loadtxt(working_directory + "/" + directory + '/parameters.txt', delimiter=',')

        #[alpha, beta, gamma_0, dist,  mu, nu, sigma]
        parameter = params[3]
        beta = 0.004811649356064012
        savedirectory = working_directory + "/" + directory + '/'
        Z_r = Z_r / np.sqrt(beta)
        Z_i = Z_i / np.sqrt(beta)

        Z = Z_r + 1j * Z_i
        Z_conj = Z_r - 1j * Z_i
        Z_modulo = np.absolute(Z)
        center_L = - 0.5 * params[3]
        center_R = 0.5 * params[3]
        J_L = np.argmin(np.abs(X - center_L))
        J_R = np.argmin(np.abs(X - center_R))
        L = X[-1] - X[0]

        Nx = len(X)
        Nt = len(T)

        ti, tf = int(0.68 * Nt), Nt
        Z_r = Z_r[ti:tf, :]
        Z_i = Z_i[ti:tf, :]
        Z_modulo = Z_modulo[ti:tf, :]
        T = T[ti:tf] - T[ti]
        dt = T[1] - T[0]

        Nx = len(X)
        Nt = len(T)

        ###      CORRELATION AND PERIOD     ###

        Z_mod_L = Z_modulo[:, 0:int(Nx / 2)]
        Z_mod_R = Z_modulo[:, int(Nx / 2) + 1:]
        X_L = X[0:int(Nx / 2)]
        X_R = X[int(Nx / 2) + 1:]

        N_L = integrate.simpson(Z_mod_L ** 2, X_L) / 6
        N_R = integrate.simpson(Z_mod_R ** 2, X_R) / 6

        CCF = np.correlate(N_L, N_R, "full")
        tau = np.arange(-T[-1], T[-1], dt)
        Ntau = len(tau)
        dtau = tau[1] - tau[0]

        CCF_max, tau_max, CCF_I = max_finder(CCF[:-1], tau, Ntau, dtau)
        ###
        fig, ((ax2), (ax1)) = plt.subplots(nrows=1, ncols=2)
        ###
        ax1.plot(tau, CCF[:-1], color="k", zorder=1)
        ax1.scatter(tau_max, CCF_max, color="r", zorder=10)
        ax1.vlines(tau_max, 0, np.amax(CCF_max) * 1.1, linestyles="--", color="r", alpha=0.3)
        ax1.grid(alpha=0.2)
        ax1.tick_params(axis='both', which='major', labelsize=15)
        ax1.set_xlabel("$t'$", fontsize=20)
        ax1.set_ylabel("$S(t')$", fontsize=20)
        ax1.set_xlim(tau[0], tau[-1])
        ax1.set_ylim(0, np.amax(CCF_max) * 1.1)
        ###
        ax2.plot(T, N_L, color="r", label="$N_L$")
        ax2.plot(T, N_R, color="b", label="$N_R$")
        ax2.set_xlabel("$t'$", fontsize=20)
        ax2.set_ylabel("$N_i$", fontsize=20)
        ax2.set_xlim(T[0], T[-1])
        ax2.set_ylim(0, np.amax(N_R) * 1.1)
        ax2.tick_params(axis='both', labelsize=15)
        ax2.legend(fontsize=15, loc="upper right")
        ax2.grid(alpha=0.2)
        plt.tight_layout()
        plt.show()
        plt.close()

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
        print(np.diff(tau_R))
        print(period)
        print(period_std)
        t_dist = np.abs(period - T)
        I_period = np.argmin(t_dist)
        I_initial = np.argmax(N_L[:I_period])
        t_initial = T[I_initial]
        t_actual = [t_initial]
        while t_actual[-1] < T[-1]:
            t_actual.append(t_actual[-1] + period)
        I = []
        for i in range(len(t_actual)):
            dist = np.abs(t_actual[i] - T)
            i_actual = np.argmin(dist)
            I.append(i_actual)

        ###      RABI OSCILLATIONS TIMESERIES     ###

        Z_timeseries_R = Z_modulo[:, J_R]
        Z_timeseries_L = Z_modulo[:, J_L]

        ###      SPATIAL PROFILES     ###
        I_final = np.argmin(np.abs(T - t_actual[1]))
        I_inicial = 0

        I_R_max = np.argmax(N_R[I_inicial:I_inicial + I_period])
        I_L_max = np.argmax(N_L[I_inicial:I_inicial + I_period])

        NR_max, TR_max, IR_max = max_finder(N_R[I_inicial:I_inicial + I_period], T[I_inicial:I_inicial + I_period], len(T[I_inicial:I_inicial + I_period]), dt)
        NL_max, TL_max, IL_max = max_finder(N_L[I_inicial:I_inicial + I_period], T[I_inicial:I_inicial + I_period], len(T[I_inicial:I_inicial + I_period]), dt)

        #plt.plot(T[I_inicial:I_inicial + I_period], N_R[I_inicial:I_inicial + I_period])
        #plt.scatter(TR_max, NR_max)
        #plt.show()
        if m == 1:
            for j in range(len(IR_max)):
                IR = IR_max[np.argmax(NR_max)]
            for j in range(len(IL_max)):
                IL = IL_max[np.argmax(NL_max)]
            profiles_R.append(Z_modulo[I_inicial + IR, :] + m * 0.5)
            profiles_L.append(Z_modulo[I_inicial + IL, :] + m * 0.5)
        else:
            dist_ZR = []
            dist_ZL = []
            for j in range(len(IR_max)):
                dist_ZR.append(np.std(profiles_R[-1] - (Z_modulo[I_inicial + IR_max[j], :] + m * 0.5)))
            for j in range(len(IL_max)):
                dist_ZL.append(np.std(profiles_L[-1] - (Z_modulo[I_inicial + IL_max[j], :] + m * 0.5)))
            print(dist_ZR)
            IR = IR_max[np.argmin(dist_ZR)]
            IL = IL_max[np.argmin(dist_ZL)]
            profiles_R.append(Z_modulo[I_inicial + IR, :] + m * 0.5)
            profiles_L.append(Z_modulo[I_inicial + IL, :] + m * 0.5)
        X_domain.append(X)

        ###      RABI OSCILLATIONS INTEGRAL LEFT/RIGHT      ###

        NL_period = N_L[I[0]:I[-1]]
        NR_period = N_R[I[0]:I[-1]]
        power_L = (1/(t_actual[-1] - t_actual[0])) * integrate.simpson((NL_period - np.mean(NL_period)) ** 2, T[I[0]:I[-1]])
        power_R = (1/(t_actual[-1] - t_actual[0])) * integrate.simpson((NR_period - np.mean(NR_period)) ** 2, T[I[0]:I[-1]])

        data = [parameter, period, period_std, power_L, power_R]
        DATA.append(data)


        ###      PLOTS      ###

        t_max = 400
        matplotlib.rc('xtick', labelsize=15)
        matplotlib.rc('ytick', labelsize=15)

        fig, (ax1, ax2, ax3) = plt.subplots(3)
        ax1.plot(T - (T[-1] - t_max), Z_modulo[:, J_R], label="$|A(x_R, t)|$", c="b", linewidth=2)
        ax1.plot(T - (T[-1] - t_max), Z_modulo[:, J_L], label="$|A(x_L, t)|$", c="r", linewidth=2)
        ax1.set_xlabel("$t/T$", fontsize=25)
        ax1.set_xlim(0, t_max)
        ax1.set_ylim(0, 5.1)
        #ax1.legend(loc="upper right", fontsize=15)
        ax1.grid(alpha=0.3)

        ax2.plot(T - (T[-1] - t_max), N_R, label="$N_{R}$", c="b", linewidth=2)
        ax2.plot(T - (T[-1] - t_max), N_L, label="$N_{L}$", c="r", linewidth=2)
        ax2.set_xlabel("$t/T$", fontsize=25)
        ax2.set_xlim(0, t_max)
        ax2.legend(loc="upper right", fontsize=15)
        ax2.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(savedirectory + '/time_series_dist.png', dpi=300)
        plt.close()

        pcm = plt.pcolormesh(X, T, Z_modulo, cmap=parula_map, shading='auto')
        cbar = plt.colorbar(pcm, shrink=1)
        cbar.set_label('$|A(x,t)|$', rotation=0, size=20, labelpad=-27, y=1.1)
        plt.xlim([X[0], X[-1]])
        plt.ylim([0, t_max])
        plt.xlabel('$x$', size='20')
        plt.ylabel('$t$', size='20')
        plt.xticks(fontsize="15")
        plt.yticks(fontsize="15")
        plt.xlim(-100, 100)
        plt.grid(linestyle='--', alpha=0.5)
        plt.savefig(savedirectory + '/module_spacetime_dist.png', dpi=300)
        plt.close()

        pcm = plt.pcolormesh(X, T, Z_modulo, cmap=parula_map, shading='auto')
        cbar = plt.colorbar(pcm, shrink=1)
        cbar.set_label('$|A(x,t)|$', rotation=0, size=20, labelpad=-27, y=1.1)
        plt.xlim([X[0], X[-1]])
        plt.xlabel('$x$', size='20')
        plt.ylabel('$t$', size='20')
        plt.xticks(fontsize="15")
        plt.yticks(fontsize="15")
        plt.xlim(-50, 50)
        plt.ylim([t_actual[0], t_actual[1]])
        plt.grid(linestyle='--', alpha=0.5)
        plt.savefig(savedirectory + '/oneperiod_spacetime_dist.png', dpi=300)
        plt.close()
    datafile = "C:/mnustes_science/simulation_data/FD/PDNLS_extended_PT/extras/dimensional/analysis"
    #fig, (ax1, ax2) = plt.subplots(2)
    #for i in range(len(directories)):
    #    ax1.plot(X_domain[i], profiles_R[i], c="r", zorder=2)
    #    ax2.plot(X_domain[i], profiles_L[i], c="b", zorder=0)
    #ax1.set_xlabel("$x\ \\textrm{(mm)}$", fontsize=20)
    #ax1.set_ylabel("$|A|\ \\textrm{(A.U.)}$", fontsize=20)
    #ax1.set_xlim([-100, 100])
    #ax2.set_xlabel("$x\ \\textrm{(mm)}$", fontsize=20)
    #ax2.set_ylabel("$|A|\ \\textrm{(A.U.)}$", fontsize=20)
    #ax2.set_xlim([-100, 100])
    #plt.tight_layout()
    #plt.savefig(datafile + '/profiles_joydivision_dist.png', dpi=300)
    #plt.close()
    DATA = np.array(DATA)
    if not os.path.exists(datafile):
        os.makedirs(datafile)
    np.savetxt(datafile + '/data_dist.txt', DATA, delimiter=',')