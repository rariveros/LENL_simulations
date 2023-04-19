from functions import *
from back_process import *
from time_integrators import *
from directories_lyap import *

if __name__ == '__main__':

    # Definiendo parámetros
    project_name = main_directory + directory_branching + "/ham_test_01"
    disc = 'E:/'
    route = 'mnustes_science/simulation_data/FD'
    eq = 'PDNLS'
    save_rate = 10
    plots = "no"
    [tmin, tmax, dt] = [0, 30000, 0.02]
    [xmin, xmax, dx] = [-70, 70, 0.2]
    d_nu = 0.0
    d_gammas = 0.002
    nus = [0.4]
    sigmas = [10] #np.arange(0.255, 0.259 + d_gammas, d_gammas)
    for nu_i in nus:
        for sigma in sigmas:
            alpha = 1
            beta = 1
            gamma_0 = 0.5
            mu_0 = 0.1
            nu = nu_i
            sigma = sigma
            [alpha_str, beta_str, mu_str, nu_str, sigma_str, gamma_str] = pdnlS_str_parameters([alpha, beta, mu_0, nu, sigma, gamma_0], 0)

            # Definición de la grilla

            t_grid = np.arange(tmin, tmax + dt, dt)
            x_grid = np.arange(xmin, xmax, dx)
            T = tmax
            Nt = t_grid.shape[0]
            Nx = x_grid.shape[0]

            # Initial Conditions Pattern
            U_1_init = 0.01 * np.random.rand(Nx)
            U_2_init = 0.01 * np.random.rand(Nx)

            # Empaquetamiento de parametros, campos y derivadas para integración
            L = xmax - xmin
            D1 = sparse_D(Nx, dx)
            D2 = sparse_DD(Nx, dx)
            operators = np.array([D2, D1])
            fields_init = [U_1_init, U_2_init]
            grids = [t_grid, x_grid, 0]
            gamma_real = gamma_0 * np.exp(- x_grid ** 2 / (2 * sigma ** 2))
            gamma_img = 0
            gamma = [gamma_real, 0]
            mu = mu_0 * np.ones(Nx)
            mu[0:10] = 10
            mu[-10:-1] = 10

            parameters = [alpha, beta, gamma, mu, nu]

            # Midiendo tiempo inicial
            now = datetime.datetime.now()
            print('Hora de Inicio: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
            time_init = time.time()

            final_fields, fields_history, time_grid = RK4_FD_ham(eq, fields_init, parameters, grids, dt, Nt, operators, save_rate)

            now = datetime.datetime.now()
            print('Hora de Término: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
            time_fin = time.time()
            print(str(time_fin - time_init) + ' seg')

            # Reobteniendo campos
            hamiltonian = np.array(fields_history)#[:, 0]

            # Guardando datos
            file = disc + route + project_name
            subfile = pdnlS_name([alpha, beta, mu_0, nu, sigma, gamma_0], "ABMNSG")
            parameters_np = np.array([alpha, beta, gamma_0, mu_0, nu])
            if not os.path.exists(file + subfile):
                os.makedirs(file + subfile)
            np.savetxt(file + subfile + '/hamiltonian.txt', hamiltonian, delimiter=',')
            np.savetxt(file + subfile + '/T.txt', time_grid, delimiter=',')
            plt.plot(time_grid, hamiltonian)
            plt.show()
            plt.close()