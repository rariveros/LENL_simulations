from time_integrators import *

if __name__ == '__main__':

    # Definiendo parámetros
    project_name = '/PT_dimer'
    disc = 'C:/'                                        # DISCO DE TRABAJO
    route = 'mnustes_science/simulation_data/FD'        # CARPETA DE TRABAJO
    eq = 'PDNLS'                                        # ECUACION
    t_rate = 1000                                      # CADA CUANTAS ITERACIONES GUARDA
    dt = 0.01
    T = 2000
    dx = 0.25#en milimetros
    ies = np.arange(0.01, 0.32, 0.02)
    jes = [0, np.pi] #np.arange(150, 50, -4)
    [tmin, tmax, dt] = [0, T, dt]
    [xmin, xmax, dx] = [-70, 70, dx]
    t_grid = np.arange(tmin, tmax + dt, dt)
    x_grid = np.arange(xmin, xmax, dx)
    T = tmax
    Nt = t_grid.shape[0]
    Nx = x_grid.shape[0]

    print("N° of simulations: " + str(len(ies) * len(jes)))
    for i in ies:
        U_1_init = 0.01 * np.random.rand(Nx)
        U_2_init = 0.01 * np.random.rand(Nx)
        for j in jes:
            alpha = 1 #6.524  #5.721
            beta = 1
            nu = i #0.014 #0.04812#0.0327449 #0.0052
            mu = 0.1
            gamma_0 = 0.28
            sigma = 3
            [alpha_str, beta_str, mu_str, nu_str, sigma_str, gamma_str] = pdnlS_str_parameters([alpha, beta, mu, nu, sigma, gamma_0], 0)

            gamma_str = str(int(gamma_0 * 1000) * 0.001)
            nu_str = str(int(nu * 1000) * 0.001)
            mu_str = str(int(mu * 1000) * 0.001)
            print('gamma = ' + gamma_str)
            print('nu = ' + nu_str)
            print('mu = ' + mu_str)

            # Empaquetamiento de parametros, campos y derivadas para integración
            L = xmax - xmin
            D2 = sparse_DD_neumann(Nx, dx)
            operators = np.array([D2])
            fields_init = [U_1_init, U_2_init]
            grids = [t_grid, x_grid, 0]

            phi = j
            gamma_real = gamma_0 * np.cos(phi) * np.exp(- (x_grid) ** 2 / (2 * sigma ** 2))
            gamma_img = gamma_0 * np.sin(phi) * np.exp(- (x_grid) ** 2 / (2 * sigma ** 2))
            gamma = [gamma_real, gamma_img]
            parameters = [alpha, beta, gamma, mu, nu]

            # Midiendo tiempo inicial
            now = datetime.datetime.now()
            print('Hora de Inicio: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
            time_init = time.time()

            final_fields, fields_history, time_grid = RK4_FD(eq, fields_init, parameters, grids, dt, Nt, operators, t_rate) #INTEGRACION EN EL TIEMPO

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

            # Guardando datos
            file = disc + route + project_name
            subfile = pdnlS_name([alpha, beta, mu, nu, sigma, gamma_0], "ABMNSG")
            parameters_np = np.array([alpha, beta, gamma_0, mu, nu, sigma, phi])
            if phi == 0:
                phi_name = "00"
            elif phi == np.pi:
                phi_name = "01"

            if not os.path.exists(file + subfile + '/phi' + phi_name):
                os.makedirs(file + subfile + '/phi' + phi_name)
            np.savetxt(file + subfile + '/phi' + phi_name + '/field_real_0.txt', U1_light[-1, :], delimiter=',')
            np.savetxt(file + subfile + '/phi' + phi_name + '/field_img_0.txt', U2_light[-1, :], delimiter=',')
            np.savetxt(file + subfile + '/phi' + phi_name + '/parameters.txt', parameters_np, delimiter=',')
            np.savetxt(file + subfile + '/phi' + phi_name + '/X.txt', x_grid, delimiter=',')

            # Guardando gráficos
            pcm = plt.pcolormesh(x_grid, t_light, modulo_light_1, cmap=parula_map, shading='auto')
            cbar = plt.colorbar(pcm, shrink=1)
            cbar.set_label('$|A|$', rotation=0, size=20, labelpad=-27, y=1.1)
            plt.xlim([x_grid[0], x_grid[-1]])
            plt.xlabel('$x$', size='20')
            plt.ylabel('$t$', size='20')
            plt.grid(linestyle='--', alpha=0.5)
            plt.title('$\gamma_0 = ' + gamma_str + '\ \\alpha = ' + alpha_str  + '\ \\textrm{mm}^{2}' + '\ \\beta = ' + beta_str + '\ \\nu = ' + nu_str + '$', size='12')
            plt.savefig(file + subfile + '/phi' + phi_name + '/module_spacetime.png', dpi=200)
            plt.close()

            pcm = plt.pcolormesh(x_grid, t_light, arg_light_1, cmap=parula_map, shading='auto')
            cbar = plt.colorbar(pcm, shrink=1)
            cbar.set_label('$\\textrm{arg}(A)$', rotation=0, size=20, labelpad=-20, y=1.1)
            plt.xlim([x_grid[0], x_grid[-1]])
            plt.xlabel('$x$', size='20')
            plt.ylabel('$t$', size='20')
            plt.grid(linestyle='--', alpha=0.5)
            plt.title(
                '$\gamma_0 = ' + gamma_str + '\ \\alpha = ' + alpha_str + '\ \\textrm{mm}^{2}' + '\ \\beta = ' + beta_str + '\ \\nu = ' + nu_str + '$',
                size='12')
            plt.savefig(file + subfile + '/phi' + phi_name + '/arg_spacetime.png', dpi=200)
            plt.close()

            pcm = plt.pcolormesh(x_grid, t_light, U1_light, cmap=parula_map, vmin=-np.amax(U1_light),
                                 vmax=np.amax(U1_light), shading='auto')
            cbar = plt.colorbar(pcm, shrink=1)
            cbar.set_label('$A_R(x, t)$', rotation=0, size=20, labelpad=-27, y=1.1)
            plt.xlim([x_grid[0], x_grid[-1]])
            plt.xlabel('$x$', size='20')
            plt.ylabel('$t$', size='20')
            plt.grid(linestyle='--', alpha=0.5)
            plt.title(
                '$\gamma_0 = ' + gamma_str + '\ \\alpha = ' + alpha_str + '\ \\textrm{mm}^{2}' + '\ \\beta = ' + beta_str + '\ \\nu = ' + nu_str + '$',
                size='12')
            plt.savefig(file + subfile + '/phi' + phi_name + '/real_spacetime.png', dpi=200)
            plt.close()

            pcm = plt.pcolormesh(x_grid, t_light, U2_light, cmap=parula_map, shading='auto')
            cbar = plt.colorbar(pcm, shrink=1)
            cbar.set_label('$A_I(x, t)$', rotation=0, size=20, labelpad=-27, y=1.1)
            plt.xlim([x_grid[0], x_grid[-1]])
            plt.xlabel('$x$', size='20')
            plt.ylabel('$t$', size='20')
            plt.grid(linestyle='--', alpha=0.5)
            plt.title(
                '$\gamma_0 = ' + gamma_str + '\ \\alpha = ' + alpha_str + '\ \\textrm{mm}^{2}' + '\ \\beta = ' + beta_str + '\ \\nu = ' + nu_str + '$',
                size='12')
            plt.savefig(file + subfile + '/phi' + phi_name + '/img_spacetime.png', dpi=200)
            plt.close()