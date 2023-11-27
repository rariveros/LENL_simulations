from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':

    # Definiendo parámetros
    project_name = '/PDNLS_scale_test'
    disc = 'C:/'
    route = 'mnustes_science/simulation_data/FD'
    eq = 'PDNLS'
    t_rate = 10
    gammas = [0.137]
    sigmas = [46]
    #directory = "D:/mnustes_science/simulation_data/FD/PDNLS_oscillatory_02/alpha=5.721/beta=1.000/mu=0.100/nu=0.027/sigma=30.000/gamma=0.190"

    #Z_r_0 = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    #Z_i_0 = np.loadtxt(directory + '/field_img.txt', delimiter=',')

    for sigma_i in sigmas:
        for gamma_i in gammas:
            alpha = 2 * 6.524
            beta = 1
            gamma_0 = gamma_i
            mu_0 = 0.1
            nu = 0.15771
            sigma = sigma_i

            # Definición de la grilla
            [tmin, tmax, dt] = [0, 5000, 0.05]
            [xmin, xmax, dx] = [-200, 200, 1]
            t_grid = np.arange(tmin, tmax + dt, dt)
            x_grid = np.arange(xmin, xmax, dx)
            T = tmax
            Nt = t_grid.shape[0]
            Nx = x_grid.shape[0]

            #np.savetxt(directory + '/field_real_init.txt', Z_r_0[-1, :], delimiter=',')
            #np.savetxt(directory + '/field_imag_init.txt', Z_i_0[-1, :], delimiter=',')


            # Initial Conditions Pattern
            #[s, c, d, phi] = special.ellipj(x_grid, 1)
            U_1_init = 0.01 * np.random.rand(Nx)#Z_r_0[-1, :]#
            U_2_init = 0.01 * np.random.rand(Nx)#Z_i_0[-1, :]#Z_i_0[-1, :]#

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

            # Guardando datos
            file = disc + route + project_name
            subfile = pdnlS_name([alpha, beta, mu_0, nu, sigma, gamma_0], "ABMNSG")
            parameters_np = np.array([alpha, beta, gamma_0, mu_0, nu, sigma])
            if not os.path.exists(file + subfile):
                os.makedirs(file + subfile)

            np.savetxt(file + subfile + '/field_real.txt', U1, delimiter=',')
            np.savetxt(file + subfile + '/field_img.txt', U2, delimiter=',')
            np.savetxt(file + subfile + '/parameters.txt', parameters_np, delimiter=',')
            np.savetxt(file + subfile + '/X.txt', x_grid, delimiter=',')
            np.savetxt(file + subfile + '/T.txt', t_grid, delimiter=',')

            pcm = plt.pcolormesh(x_grid, t_light[::10], modulo_light[::10, :]/np.sqrt(0.004811649356064012), cmap=parula_map, shading='auto')
            cbar = plt.colorbar(pcm, shrink=1)
            cbar.set_label('$|A|$', rotation=0, size=20, labelpad=-27, y=1.1)
            plt.xlim([x_grid[0], x_grid[-1]])
            plt.xlabel('$x$', size='20')
            plt.ylabel('$t$', size='20')
            plt.grid(linestyle='--', alpha=0.5)
            plt.title(
                '$\gamma = ' + str(gamma_str) + '$' + '    $\\nu = ' + str(nu_str) + '$' + '    $\mu = ' + str(mu_str) + '$',
                size='20')
            plt.savefig(file + subfile + '/module_spacetime.png', dpi=300)
            plt.close()

            pcm = plt.pcolormesh(x_grid, t_light[::10], arg_light[::10, :], cmap=parula_map, shading='auto')
            cbar = plt.colorbar(pcm, shrink=1)
            cbar.set_label('$\\textrm{arg}(A)$', rotation=0, size=20, labelpad=-20, y=1.1)
            plt.xlim([x_grid[0], x_grid[-1]])
            plt.xlabel('$x$', size='20')
            plt.ylabel('$t$', size='20')
            plt.grid(linestyle='--', alpha=0.5)
            plt.title(
                '$\gamma = ' + str(gamma_str) + '$' + '    $\\nu = ' + str(nu_str) + '$' + '    $\mu = ' + str(mu_str) + '$',
                size='20')
            plt.savefig(file + subfile + '/arg_spacetime.png', dpi=300)
            plt.close()

            pcm = plt.pcolormesh(x_grid, t_light[::10], U1_light[::10, :], cmap=parula_map,vmin=-np.amax(U1_light), vmax=np.amax(U1_light), shading='auto')
            cbar = plt.colorbar(pcm, shrink=1)
            cbar.set_label('$A_R(x, t)$', rotation=0, size=20, labelpad=-27, y=1.1)
            plt.xlim([x_grid[0], x_grid[-1]])
            plt.xlabel('$x$', size='20')
            plt.ylabel('$t$', size='20')
            plt.grid(linestyle='--', alpha=0.5)
            plt.title(
                '$\gamma = ' + str(gamma_str) + '$' + '    $\\nu = ' + str(nu_str) + '$' + '    $\mu = ' + str(mu_str) + '$',
                size='20')
            plt.savefig(file + subfile + '/real_spacetime.png', dpi=300)
            plt.close()

            pcm = plt.pcolormesh(x_grid, t_light[::10], U2_light[::10, :], cmap=parula_map, shading='auto')
            cbar = plt.colorbar(pcm, shrink=1)
            cbar.set_label('$A_I(x, t)$', rotation=0, size=20, labelpad=-27, y=1.1)
            plt.xlim([x_grid[0], x_grid[-1]])
            plt.xlabel('$x$', size='20')
            plt.ylabel('$t$', size='20')
            plt.grid(linestyle='--', alpha=0.5)
            plt.title(
                '$\gamma = ' + str(gamma_str) + '$' + '    $\\nu = ' + str(nu_str) + '$' + '    $\mu = ' + str(mu_str) + '$',
                size='20')
            plt.savefig(file + subfile + '/img_spacetime.png', dpi=300)
            plt.close()

            plt.plot(x_grid, modulo_light[-1, :]/np.sqrt(0.004811649356064012), label="$|A(x,t)|$")
            plt.xlabel('$x$', size='25')
            plt.ylabel('$|A(x)|$', size='25')
            plt.legend(fontsize=18)
            plt.xlim([-180, 180])
            plt.grid(alpha=0.3)
            plt.tight_layout()
            plt.savefig(file + subfile + '/final_profile.png', dpi=200)
            plt.close()