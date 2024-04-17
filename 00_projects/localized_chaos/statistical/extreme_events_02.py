from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':

    # Definiendo parámetros
    project_name = '/localized_chaos/EE_02'
    disc = 'C:/'
    route = 'mnustes_science/simulation_data/FD'
    eq = 'PDNLS'
    t_rate = 1
    [tmin, tmax, dt] = [0, 4000, 0.02]
    [xmin, xmax, dx] = [-100, 100, 0.2]
    t_grid = np.arange(tmin, tmax + dt, dt)
    x_grid = np.arange(xmin, xmax, dx)
    T = tmax
    Nt = t_grid.shape[0]
    Nx = x_grid.shape[0]
    D_kys = []
    param = []
    U_1_init = 0.01 * np.random.rand(Nx)  # Z_r_0[-1, :]#
    U_2_init = 0.01 * np.random.rand(Nx)  # Z_i_0[-1, :]#Z_i_0[-1, :]#
    ies = np.arange(0.275, 0.525, 0.025)#
    jotas = [10]
    m = 100.18
    i = 0
    Z_MAX = []
    PDF_EE = []
    dH_EE = []
    P = []
    file = disc + route + project_name
    n = 0
    if not os.path.exists(file):
        os.makedirs(file)
    for gamma_i in ies:
        for sigma_i in jotas:
            n = n + 1
            print("###### N° of simulations: " + str(n) + "/" + str(len(ies)*len(jotas)) + " ######")
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
            U_1_init = U1[-1, :]
            U_2_init = U2[-1, :]
            [tmin, tmax, dt] = [0, 1500, 0.02]
            T = tmax
            t_grid = np.arange(tmin, tmax + dt, dt)
            Nt = t_grid.shape[0]
            if gamma_i != ies[0]:
                Z_video = modulo[-50000::, :]
                T_video = time_grid[-50000::]
                Nt_video = len(T_video)
                Z_maxs = []
                X_maxs = []
                for j in range(Nx):
                    if x_grid[j] > -40 and x_grid[j] < 40:
                        Z_max, X_max, I_max = max_finder(Z_video[:, j], T_video, Nt_video, dt)
                        Z_maxs.extend(Z_max)
                        X_maxs.extend(X_max)
                Z_maxs = np.array(Z_maxs)
                X_maxs = np.array(X_maxs)
                #plt.plot(x_grid, Z_video[-1, :])
                #plt.vlines([-45, 45], 0, 1.2 * np.amax(Z_video[-1, :]), colors="r", linestyles="--")
                #plt.ylim([0, 1.1 * np.amax(Z_video[-1, :])])
                #plt.show()
                hist, bins = np.histogram(Z_maxs, 200, range=(0., 4))
                hist = hist / np.sum(hist)
                Q1 = np.quantile(Z_maxs, .33333)
                Q2 = np.quantile(Z_maxs, .66667)
                Z_maxs_q2 = Z_maxs[Z_maxs > Q2]
                H_3 = np.mean(Z_maxs_q2)
                H_mean = np.mean(Z_maxs)
                #EE = Z_maxs[Z_maxs > 2 * H_3]
                EE = Z_maxs[Z_maxs > 3 * H_mean]
                p_EE = len(EE) / len(Z_maxs)
                P.append(p_EE)
                print("N° Extreme Events = " + str(len(EE)))
                dH = []
                for i in range(1, len(hist) + 1):
                    dH.append((bins[i] + bins[i - 1]) / 2)
                #plt.vlines([Q1, Q2], 0, 1.2 * np.amax(hist), colors="r", linestyles="--")
                #plt.vlines(3 * H_mean, -12, 0, colors="k", linestyles="--")
                #plt.ylim([0, 1.1 * np.amax(hist)])
                #plt.plot(dH, np.log(hist))
                #plt.show()
                #plt.close()
                Z_MAX.append(Z_maxs)
                PDF_EE.append(np.log(hist))
                dH_EE.append(dH)
    #np.savetxt(file + '/Z_maxs.txt', Z_MAX, delimiter=',')
    np.savetxt(file + '/PDF_EEs.txt', PDF_EE, delimiter=',')
    np.savetxt(file + '/dH_EEs.txt', dH_EE, delimiter=',')
    np.savetxt(file + '/p_EE.txt', P, delimiter=',')
    np.savetxt(file + '/gammas.txt', ies, delimiter=',')
    plt.scatter(ies, P)
    plt.show()
    plt.close()


