from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':

    # Definiendo parámetros
    project_name = '/PDNLS_oscillatory_02/diff'
    disc = 'D:/'
    route = 'mnustes_science/simulation_data/FD'
    eq = 'PDNLS'
    t_rate = 10
    gammas = [0.19]
    sigmas = [30]
    directory = "D:/mnustes_science/simulation_data/FD/PDNLS_oscillatory_02/alpha=5.721/beta=1.000/mu=0.100/nu=0.027/sigma=30.000/gamma=0.190"

    Z_r_0 = np.loadtxt(directory + '/field_real_init.txt', delimiter=',')
    Z_i_0 = np.loadtxt(directory + '/field_imag_init.txt', delimiter=',')

    for sigma_i in sigmas:
        for gamma_i in gammas:
            alpha = 5.721
            beta = 1
            gamma_0 = gamma_i
            mu_0 = 0.1
            nu = 0.0274
            sigma = sigma_i

            # Definición de la grilla
            [tmin, tmax, dt] = [0, 500, 0.01]
            [xmin, xmax, dx] = [-200, 200, 0.5]
            t_grid = np.arange(tmin, tmax + dt, dt)
            x_grid = np.arange(xmin, xmax, dx)
            T = tmax
            Nt = t_grid.shape[0]
            Nx = x_grid.shape[0]

            # Initial Conditions Pattern
            #[s, c, d, phi] = special.ellipj(x_grid, 1)
            U_1_init = Z_r_0 #+ 0.01 * np.random.rand(Nx)#
            U_2_init = Z_i_0 #+ 0.01 * np.random.rand(Nx)#Z_i_0[-1, :]
            #U_1_init[int(Nx / 2)] = U_1_init[int(Nx / 2)] + 0.01 * 2 * (np.random.rand(1) - 0.5)
            #U_2_init[int(Nx / 2)] = U_2_init[int(Nx / 2)] + 0.01 * 2 * (np.random.rand(1) - 0.5)
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
            mu[0:10] = 10
            mu[-10:-1] = 10

            parameters = [alpha, beta, gamma, mu_0, nu]
            [alpha_str, beta_str, mu_str, nu_str, sigma_str, gamma_str] = pdnlS_str_parameters([alpha, beta, mu_0, nu, sigma, gamma_0], 0)

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

            lightness = 3
            U1_light = U1[0:-1:lightness]
            U2_light = U2[0:-1:lightness]
            U_complex_light = U_complex[0:-1:lightness]
            t_light = time_grid[0:-1:lightness]
            modulo_light = modulo[0:-1:lightness]
            arg_light = arg[0:-1:lightness]

            dU1 = np.abs(U1 - Z_r_0)
            dU2 = np.abs(U2 - Z_i_0)
            dU_complex = dU1 + 1j * dU2
            dU_module = np.absolute(dU_complex)
            darg = np.angle(dU_complex)
            darg = (2 * np.pi + darg) * (darg < 0) + darg * (darg > 0)




            # Guardando datos
            file = disc + route + project_name
            subfile = pdnlS_name([alpha, beta, mu_0, nu, sigma, gamma_0], "ABMNSG")
            parameters_np = np.array([alpha, beta, gamma_0, mu_0, nu, sigma])
            if not os.path.exists(file + subfile):
                os.makedirs(file + subfile)

            pcm = plt.pcolormesh(x_grid, t_grid[::10], dU_module[::10, :], cmap=parula_map, shading='auto')
            cbar = plt.colorbar(pcm, shrink=1)
            cbar.set_label('$|\delta A|$', rotation=0, size=20, labelpad=-27, y=1.1)
            plt.xlim([x_grid[0], x_grid[-1]])
            plt.xlabel('$x$', size='20')
            plt.ylabel('$t$', size='20')
            plt.grid(linestyle='--', alpha=0.5)
            plt.title(
                '$\gamma = ' + str(gamma_str) + '$' + '    $\\nu = ' + str(nu_str) + '$' + '    $\mu = ' + str(mu_str) + '$',
                size='20')
            plt.savefig(file + subfile + '/module_spacetime.png', dpi=300)
            plt.close()