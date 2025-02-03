import matplotlib.pyplot as plt
from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':

    # Definiendo parámetros
    project_name = '/soliton_reduced'
    disc = 'D:/'
    route = 'mnustes_science/simulation_data/FD'
    eq = 'reduced_soliton'
    t_rate = 10

    [alpha, beta, mu, sigma, gamma] = [6.524, 1, 0.075, 15, 0.18]
    alpha_str = f"{alpha:.{3}f}"
    beta_str = f"{beta:.{3}f}"
    mu_str = f"{mu:.{3}f}"
    sigma_str = f"{sigma:.{2}f}"
    gamma_str = f"{gamma:.{3}f}"

    working_directory = 'D:/mnustes_science/simulation_data/FD/soliton_reduced/alpha=' + alpha_str + '/beta=' + beta_str + '/mu=' + mu_str
    directories = [name for name in os.listdir(working_directory) if os.path.isdir(os.path.join(working_directory, name))]
    for directory in directories:
        dir = working_directory + "/" + directory
        D1 = np.loadtxt(dir + '/gamma=' + gamma_str + '/sigma=' + sigma_str + '/D1.txt', delimiter=',')
        D2 = np.loadtxt(dir + '/gamma=' + gamma_str + '/sigma=' + sigma_str + '/D2.txt', delimiter=',')
        parameters = np.loadtxt(dir + '/gamma=' + gamma_str + '/sigma=' + sigma_str + '/parameters.txt', delimiter=',')

        ti = 0.

        # Definición de la grilla
        [tmin, tmax, dt] = [0, 5000, 1]
        [xmin, xmax, dx] = [0, 1, 1]
        t_grid = np.arange(tmin, tmax + dt, dt)
        x_grid = np.arange(xmin, xmax, dx)
        T = tmax
        Nt = t_grid.shape[0]
        Nx = x_grid.shape[0]

        # Initial Conditions
        U_init = np.ones(Nx)
        V_init = np.ones(Nx)
        operators = [D1, D2]

        # Empaquetamiento de parametros, campos y derivadas para integración
        fields_init = [U_init, V_init]
        grids = [t_grid, x_grid, 0]
        parameters_np = parameters #[alpha, beta, mu, nu, sigma, gamma]
        [alpha, beta, mu, nu, sigma, gamma] = parameters
        # Midiendo tiempo inicial
        now = datetime.datetime.now()
        print('Hora de Inicio: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
        time_init = time.time()

        final_fields, fields_history, time_grid = RK4_FD(eq, fields_init, parameters_np, grids, dt, Nt, operators, t_rate)

        now = datetime.datetime.now()
        print('Hora de Término: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
        time_fin = time.time()
        print(str(time_fin - time_init) + ' seg')

        # Reobteniendo campos
        U = np.array(fields_history)[:, 0]
        V = np.array(fields_history)[:, 1]

        # Guardando datos
        file = disc + route + project_name
        subfile = "/test"
        if not os.path.exists(file + subfile):
            os.makedirs(file + subfile)

        np.savetxt(file + subfile + '/U.txt', U, delimiter=',')
        np.savetxt(file + subfile + '/V.txt', V, delimiter=',')
        np.savetxt(file + subfile + '/parameters.txt', parameters_np, delimiter=',')
        np.savetxt(file + subfile + '/T.txt', t_grid, delimiter=',')

        lightness = 1
        U_light = U[0::lightness]
        V_light = V[0::lightness]
        t_light = time_grid[0::lightness]

        #fig, (ax1) = plt.subplots(1)
        #ax1.plot(t_light, U_light, c="b")
        #ax1.plot(t_light, V_light, c="r")
        #plt.xlabel('$x$', size='20')
        #plt.ylabel('$t$', size='20')
        #ax1.set_xlim([ti * t_light[-1], t_light[-1]])
        #ax1.grid(linestyle='--', alpha=0.5)
        plt.scatter(nu, U_light[-1], c="b")
        plt.scatter(nu, V_light[-1], c="r")
    plt.show()
    plt.close()





