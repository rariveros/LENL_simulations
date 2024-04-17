from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':

    # Definiendo parámetros
    project_name = '/localized_chaos/EE'
    disc = 'C:/'
    route = 'mnustes_science/simulation_data/FD'
    eq = 'PDNLS'
    t_rate = 1
    [tmin, tmax, dt] = [0, 2000, 0.02]
    [xmin, xmax, dx] = [-100, 100, 0.2]
    t_grid = np.arange(tmin, tmax + dt, dt)
    x_grid = np.arange(xmin, xmax, dx)
    T = tmax
    Nt = t_grid.shape[0]
    Nx = x_grid.shape[0]
    D_kys = []
    param = []
    N_condit = int(1 * 2 * Nx)
    x_condit = np.arange(0, N_condit)/(2*Nx)
    U_1_init = 0.01 * np.random.rand(Nx)  # Z_r_0[-1, :]#
    U_2_init = 0.01 * np.random.rand(Nx)  # Z_i_0[-1, :]#Z_i_0[-1, :]#
    ies = np.arange(0.2, 0.5, 0.025)#
    jotas = [10]
    m = 10
    i = 0
    for gamma_i in ies:
        for sigma_i in jotas:
            print("N° of simulations: " + str(len(ies)*len(jotas)))
            # Definición de parámetros
            alpha = 1
            beta = 1
            gamma_0 = gamma_i
            param.append(gamma_0)
            mu_0 = 0.1
            nu = 0.189
            sigma = sigma_i

            # Empaquetamiento de parametros, campos y derivadas para integración
            L = xmax - xmin
            D2 = sparse_DD_neumann(Nx, dx)
            operators = np.array([D2])
            fields_init = [U_1_init, U_2_init]
            grids = [t_grid, x_grid, 0]

            gamma_real = gamma_0 * np.exp(- x_grid ** 2 / (2 * sigma ** 2))
            gamma_img = 0
            gamma = [gamma_real, gamma_img]
            mu = mu_0 * np.ones(Nx)
            mu[0:25] = 10
            mu[-25:] = 10

            parameters = [alpha, beta, gamma, mu, nu]
            [alpha_str, beta_str, mu_str, nu_str, sigma_str, gamma_str] = pdnlS_str_parameters([alpha, beta, mu_0, nu, sigma, gamma_0], 0)

            # SIMULANDO Y CALCULANDO ESPECTRO
            # Midiendo tiempo inicial
            now = datetime.datetime.now()
            print('Hora de Inicio: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
            time_init = time.time()

            final_fields, fields_history, time_grid = RK4_FD(eq, fields_init, parameters, grids, dt, Nt, operators, t_rate)

            now = datetime.datetime.now()
            print('Hora de Término: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
            time_fin = time.time()
            print(str(time_fin - time_init) + ' seg')

            # Reobteniendo campos
            U1 = np.array(fields_history)[:, 0]
            U2 = np.array(fields_history)[:, 1]
            U_complex = U1 + 1j * U2
            t_grid = time_grid
            modulo = np.absolute(U_complex)
            arg = np.angle(U_complex)
            arg = (2 * np.pi + arg) * (arg < 0) + arg * (arg > 0)

            Z_video = modulo[-100000::, :]
            T_video = time_grid[-100000::]
            Nt_video = len(T_video)
            Z_maxs = []
            X_maxs = []
            for i in range(Nt_video):
                Z_max, X_max, J_max = max_finder(Z_video[i, :], x_grid, len(x_grid), dx)
                Z_maxs.extend(Z_max)
                X_maxs.extend(X_max)
            Z_maxs = np.array(Z_maxs)
            X_maxs = np.array(X_maxs)
            plt.plot(x_grid, Z_video[-1, :])
            plt.vlines([-45, 45], 0, 1.2 * np.amax(Z_video[-1, :]), colors="r", linestyles="--")
            plt.ylim([0, 1.1 * np.amax(Z_video[-1, :])])
            plt.show()
            Z_maxs = Z_maxs[(X_maxs >= -45) & (X_maxs <= 45)]
            hist, bins = np.histogram(Z_maxs, 200, range=(0, 4))
            Q2 = np.quantile(Z_maxs, .66667)
            Z_maxs_q2 = Z_maxs[Z_maxs > Q2]
            H_3 = np.mean(Z_maxs_q2)
            EE = Z_maxs[Z_maxs > 2 * H_3]
            print("N° Extreme Events = " + str(len(EE)))
            dH = []
            for i in range(1, len(hist) + 1):
                dH.append((bins[i] + bins[i - 1]) / 2)
            plt.vlines(Q2, 0, 2, colors="r", linestyles="--")
            plt.vlines(2 * H_3, 0, 2, colors="k", linestyles="--")
            plt.ylim([0, 1.5])
            plt.plot(dH, hist)
            plt.show()
            plt.close()
