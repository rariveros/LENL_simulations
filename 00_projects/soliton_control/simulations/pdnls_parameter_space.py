from functions import *
from back_process import *
from time_integrators import *
#from directories_lyap import *

if __name__ == '__main__':

    # Definiendo parámetros
    project_name = "/soliton_draft/variando_sigma"
    disc = 'C:/'
    route = 'mnustes_science/simulation_data/FD'
    eq = 'PDNLS'
    save_rate = 250
    plots = "si"
    file = disc + route + project_name
    [tmin, tmax, dt] = [0, 2000, 0.005]
    [xmin, xmax, dx] = [-120, 120, 0.5]
    t_grid = np.arange(tmin, tmax + dt, dt)
    x_grid = np.arange(xmin, xmax, dx)
    T = tmax
    Nt = t_grid.shape[0]
    Nx = x_grid.shape[0]

    d = 20
    f_i = np.arange(12.9, 13.3, 0.05)
    a = 0.180 / (f_i ** 2 * 8.401093 * 10 ** (-5))
    alpha, beta, nus, gamma_0, f_0 = fluid_pdnls_parameters(f_i, a, d)
    gamma_0 = gamma_0[0]
    gammas = [0.18]
    nus = np.arange(-0.15, -0.05, 0.005)
    sigmas = [12.5, 17.5]
    t_0 = tmax
    x_0 = -1
    X_max = []
    for sigma in sigmas:
        print("###### gamma = " + str(gamma_0) + " ######")
        for nu in nus:
            print("### nu = " + str(nu) + " ###")
            alpha = 4 * 6.524  # 5.721
            beta = 1
            mu_0 = 0.075
            #sigma = 15
            gamma = 0.18
            delta = np.sqrt(- nu + np.sqrt(gamma_0 ** 2 - mu_0 ** 2))
            [alpha_str, beta_str, mu_str, nu_str, sigma_str, gamma_str] = pdnlS_str_parameters([alpha, beta, mu_0, nu, sigma, gamma_0], 0)


            if nu == nus[0]:
                U_1_init = (np.sqrt(2) * delta / np.cosh(delta * (x_grid - x_0) / np.sqrt(alpha))) * np.real(np.exp(-1j * 0.5 * np.arcsin(mu_0 / gamma_0)))
                U_2_init = (np.sqrt(2) * delta / np.cosh(delta * (x_grid - x_0) / np.sqrt(alpha))) * np.imag(np.exp(-1j * 0.5 * np.arcsin(mu_0 / gamma_0)))
                [tmin, tmax, dt] = [0, 2000, 0.005]
                [xmin, xmax, dx] = [-120, 120, 0.5]
                t_grid = np.arange(tmin, tmax + dt, dt)
                x_grid = np.arange(xmin, xmax, dx)
                T = tmax
                Nt = t_grid.shape[0]
                Nx = x_grid.shape[0]

            # Empaquetamiento de parametros, campos y derivadas para integración
            L = xmax - xmin
            D2 = sparse_DD(Nx, dx)
            operators = np.array([D2])
            fields_init = [U_1_init, U_2_init]
            grids = [t_grid, x_grid, 0]
            gamma_real = gamma_0 * np.exp(- x_grid ** 2 / (2 * sigma ** 2))
            gamma_img = gamma_0 * np.exp(- x_grid ** 2 / (2 * sigma ** 2)) * 0
            gamma = [gamma_real, gamma_img]
            mu = mu_0 * np.ones(Nx)
            mu[0:10] = 10
            mu[-10:-1] = 10

            parameters = [alpha, beta, gamma, mu, nu]

            # Midiendo tiempo inicial
            now = datetime.datetime.now()
            print('Hora de Inicio: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
            time_init = time.time()

            final_fields, fields_history, time_grid = RK4_FD(eq, fields_init, parameters, grids, dt, Nt, operators, save_rate)

            now = datetime.datetime.now()
            print('Hora de Término: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
            time_fin = time.time()
            print(str(time_fin - time_init) + ' seg')

            # Reobteniendo campos
            U1_light = np.array(fields_history)[:, 0]
            U2_light = np.array(fields_history)[:, 1]
            U_complex = U1_light + 1j * U2_light
            t_light = time_grid

            # Definiendo variables finales
            modulo_light_1 = np.absolute(U_complex)
            arg_light_1 = np.angle(U_complex)
            arg_light_1 = (2 * np.pi + arg_light_1) * (arg_light_1 < 0) + arg_light_1 * (arg_light_1 > 0)
            analytical_signal_1 = hilbert(U1_light[-1, :])
            amplitude_envelope_1 = np.abs(analytical_signal_1)

            Z_last = modulo_light_1[-1, :]
            J_max = np.argmax(Z_last)
            delta_wind = 5
            window_L = J_max - delta_wind
            window_R = J_max + delta_wind

            def quadratic(x, a, b, c):
                return a * x ** 2 + b * x + c

            def soliton(x, a, b):
                return np.sqrt(2) * a / np.cosh(b * x )

            popt, pcov = curve_fit(quadratic, x_grid[window_L:window_R], Z_last[window_L:window_R])
            [a, b, c] = popt
            center = - b / (2 * a)

            popt, pcov = curve_fit(soliton, x_grid - center, Z_last)
            [amplitud, ancho] = popt
            [amplitud, ancho] = [amplitud, ancho * np.sqrt(alpha)]
            data = [center, delta, ancho, amplitud]

            # Guardando datos

            subfile = pdnlS_name([alpha, beta, mu_0, nu, sigma, gamma_0], "ABMNSG")
            parameters_np = np.array([alpha, beta, gamma_0, mu_0, nu, sigma])

            U_1_init = U1_light[-1, :]
            U_2_init = U2_light[-1, :]

            if not os.path.exists(file + subfile):
                os.makedirs(file + subfile)
            np.savetxt(file + subfile + '/data.txt', data, delimiter=',')
            np.savetxt(file + subfile + '/field_real.txt', U1_light, delimiter=',')
            np.savetxt(file + subfile + '/field_img.txt', U2_light, delimiter=',')
            np.savetxt(file + subfile + '/parameters.txt', parameters_np, delimiter=',')
            np.savetxt(file + subfile + '/final_envelope.txt', amplitude_envelope_1, delimiter=',')
            np.savetxt(file + subfile + '/X.txt', x_grid, delimiter=',')
            np.savetxt(file + subfile + '/T.txt', t_light, delimiter=',')

            if plots == "si":
                pcm = plt.pcolormesh(x_grid, t_light, modulo_light_1, cmap=parula_map, shading='auto')
                cbar = plt.colorbar(pcm, shrink=1)
                cbar.set_label('$|A|$', rotation=0, size=20, labelpad=-27, y=1.1)
                cbar.ax.tick_params(labelsize=15)
                plt.xlim([x_grid[0], x_grid[-1]])
                plt.xlabel('$x$', size='25')
                plt.ylabel('$t$', size='25')
                plt.xticks(fontsize=15)
                plt.yticks(fontsize=15)
                #plt.xlim([-50, 50])
                plt.grid(linestyle='--', alpha=0.2, color='k')
                plt.tight_layout()
                plt.savefig(file + subfile + '/module_spacetime.png', dpi=200)
                plt.close()

                pcm = plt.pcolormesh(x_grid, t_light, arg_light_1, cmap=parula_map, shading='auto')
                cbar = plt.colorbar(pcm, shrink=1)
                cbar.set_label('$\\textrm{arg}(A)$', rotation=0, size=20, labelpad=-20, y=1.1)
                cbar.ax.tick_params(labelsize=15)
                plt.xlim([x_grid[0], x_grid[-1]])
                plt.xlabel('$x$', size='25')
                plt.ylabel('$t$', size='25')
                plt.xticks(fontsize=15)
                plt.yticks(fontsize=15)
                plt.grid(linestyle='--', alpha=0.2, color='k')
                plt.tight_layout()
                plt.savefig(file + subfile + '/arg_spacetime.png', dpi=200)
                plt.close()

                pcm = plt.pcolormesh(x_grid, t_light, U1_light, cmap=parula_map, shading='auto')
                cbar = plt.colorbar(pcm, shrink=1)
                cbar.set_label('$A_R(x, t)$', rotation=0, size=20, labelpad=-27, y=1.1)
                cbar.ax.tick_params(labelsize=15)
                plt.xlim([x_grid[0], x_grid[-1]])
                plt.xlabel('$x$', size='25')
                plt.ylabel('$t$', size='25')
                plt.xticks(fontsize=15)
                plt.yticks(fontsize=15)
                plt.grid(linestyle='--', alpha=0.2, color='k')
                plt.tight_layout()
                plt.savefig(file + subfile + '/real_spacetime.png', dpi=200)
                plt.close()

                pcm = plt.pcolormesh(x_grid, t_light, U2_light, cmap=parula_map, shading='auto')
                cbar = plt.colorbar(pcm, shrink=1)
                cbar.set_label('$A_I(x, t)$', rotation=0, size=20, labelpad=-27, y=1.1)
                cbar.ax.tick_params(labelsize=15)
                plt.xlim([x_grid[0], x_grid[-1]])
                plt.xlabel('$x$', size='25')
                plt.ylabel('$t$', size='25')
                plt.xticks(fontsize=15)
                plt.yticks(fontsize=15)
                plt.grid(linestyle='--', alpha=0.2, color='k')
                plt.tight_layout()
                plt.savefig(file + subfile + '/img_spacetime.png', dpi=200)
                plt.close()

                plt.plot(x_grid, gamma_real, label="$\gamma_R$")
                plt.plot(x_grid, gamma_img, label="$\gamma_I$")
                plt.xlabel('$x$', size='25')
                plt.ylabel('$\gamma(x)$', size='25')
                plt.legend(fontsize=18)
                plt.savefig(file + subfile + '/forcing_01.png', dpi=200)
                plt.close()

                plt.plot(x_grid, modulo_light_1[-1, :], label="$A_f$")
                plt.plot(x_grid, modulo_light_1[0, :], label="$A_i$")
                plt.xlabel('$x$', size='25')
                plt.ylabel('$A(x)$', size='25')
                plt.legend(fontsize=18)
                plt.xlim([-50, 50])
                plt.grid(alpha=0.3)
                plt.tight_layout()
                plt.savefig(file + subfile + '/final_profile.png', dpi=200)
                plt.close()

                plt.plot(x_grid, np.unwrap(arg_light_1[-1, :], period=2 * np.pi))
                plt.xlabel('$x$', size='25')
                plt.ylabel('$\\theta(x)$', size='25')
                plt.tight_layout()
                plt.savefig(file + subfile + '/final_phase.png', dpi=200)
                plt.close()

            [tmin, tmax, dt] = [0, 500, 0.005]
            [xmin, xmax, dx] = [-120, 120, 0.5]
            t_grid = np.arange(tmin, tmax + dt, dt)
            x_grid = np.arange(xmin, xmax, dx)
            T = tmax
            Nt = t_grid.shape[0]
            Nx = x_grid.shape[0]





