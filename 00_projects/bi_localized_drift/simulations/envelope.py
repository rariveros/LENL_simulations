from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':

    # Definiendo parámetros
    project_name = '/bi_localized_drift/one_envelope'
    disc = 'D:/'                                        # DISCO DE TRABAJO
    route = 'mnustes_science/simulation_data/FD'        # CARPETA DE TRABAJO
    eq = 'amplitude_one_envelope'                         # ECUACION
    t_rate = 10                                         # CADA CUANTAS ITERACIONES GUARDA
    dt = 0.02
    T = 2500
    dx = 1 #en milimetros
    ies = np.arange(0.19, 0.20, 0.01)
    jotas = [1] #np.arange(0.25, 0.18, - 0.005)

    [tmin, tmax, dt] = [0, T, dt]
    [xmin, xmax, dx] = [-150, 150, dx]
    t_grid = np.arange(tmin, tmax + dt, dt)
    x_grid = np.arange(xmin, xmax, dx)
    T = tmax
    Nt = t_grid.shape[0]
    Nx = x_grid.shape[0]

    print("N° of simulations: " + str(len(ies) * len(jotas)))
    for i in ies:
        U_1_init = 0.01 * np.random.rand(Nx)
        for j in jotas:
            alpha = 6.524  #5.721
            beta = 1
            nu = 0.1 #0.014 #0.04812#0.0327449 #0.0052
            mu = 0.1
            gamma_0 = i
            dist = 80
            sigma = 16
            [alpha_str, beta_str, mu_str, nu_str, sigma_str, gamma_str] = pdnlS_str_parameters([alpha, beta, mu, nu, sigma, gamma_0], 0)

            gamma_str = str(int(gamma_0 * 1000) * 0.001)
            nu_str = str(int(nu * 1000) * 0.001)
            mu_str = str(int(mu * 1000) * 0.001)
            print('gamma = ' + gamma_str)
            print('nu = ' + nu_str)
            print('mu = ' + mu_str)
            print("dist = " + str(dist))

            # Empaquetamiento de parametros, campos y derivadas para integración
            L = xmax - xmin
            D1 = sparse_D_neumann(Nx, dx)
            D2 = sparse_DD_neumann(Nx, dx)
            operators = np.array([D1, D2])
            fields_init = [U_1_init]
            grids = [t_grid, x_grid, 0]
            phi = np.pi
            gamma_real = gamma_0
            gamma_img = 0 * gamma_0
            gamma = [gamma_real, gamma_img]

            parameters = [alpha, beta, gamma, mu, nu, sigma, dist]

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
            U1_light = np.array(fields_history)[:, 0]
            t_light = time_grid

            DU1 = []
            for i in range(len(t_light)):
                DU1.append(np.append(np.diff(U1_light[i, :]), 0) / dx)
            DU1 = np.array(DU1)

            K1c = np.sqrt(nu / alpha) #+ (gamma_0 - mu) * np.abs(U1_light) ** 2
            K2c = np.sqrt(nu / alpha) #+ (gamma_0 - mu) * np.abs(U2_light) ** 2

            psi_1 = U1_light * np.exp(1j * K1c * x_grid)
            psi_2a = (-1 / (2 * mu))* (3 * 1j * U1_light ** 3 - 2 * np.sqrt(nu / alpha) * DU1) * np.exp(1j * K1c * x_grid)
            psi_2b = (1 / (8 * nu)) * U1_light ** 3 * np.exp(3 * 1j * K1c * x_grid)
            psi = psi_1 + np.conjugate(psi_1) + psi_2a - np.conjugate(psi_2a) + psi_2b + np.conjugate(psi_2b)

            # Definiendo variables finales
            modulo_light_1 = np.absolute(U1_light)

            # Guardando datos
            file = disc + route + project_name
            subfile = pdnlS_bigauss_name([alpha, beta, mu, nu, sigma, gamma_0, dist], "ABMNSGD")
            parameters_np = np.array([alpha, beta, gamma_0, dist,  mu, nu, sigma])

            if not os.path.exists(file + subfile):
                os.makedirs(file + subfile)
            np.savetxt(file + subfile + '/field.txt', U1_light, delimiter=',')
            np.savetxt(file + subfile + '/parameters.txt', parameters_np, delimiter=',')
            np.savetxt(file + subfile + '/X.txt', x_grid, delimiter=',')
            np.savetxt(file + subfile + '/T.txt', t_light, delimiter=',')

            #################### Guardando Gráficos ####################
            pcm = plt.pcolormesh(x_grid, t_light, modulo_light_1, cmap=parula_map, shading='auto')
            cbar = plt.colorbar(pcm, shrink=1)
            cbar.set_label('$|A|$', rotation=0, size=20, labelpad=-27, y=1.1)
            plt.xlim([x_grid[0], x_grid[-1]])
            plt.xlabel('$x$', size='20')
            plt.ylabel('$t$', size='20')
            plt.grid(linestyle='--', alpha=0.5)
            plt.savefig(file + subfile + '/A_module.png', dpi=300)
            plt.close()

            pcm = plt.pcolormesh(x_grid, t_light, np.abs(psi), cmap=parula_map, shading='auto')
            cbar = plt.colorbar(pcm, shrink=1)
            cbar.set_label('$|\psi|$', rotation=0, size=20, labelpad=-27, y=1.1)
            plt.xlim([x_grid[0], x_grid[-1]])
            plt.xlabel('$x$', size='20')
            plt.ylabel('$t$', size='20')
            plt.grid(linestyle='--', alpha=0.5)
            plt.savefig(file + subfile + '/pdnlS_module_spacetime.png', dpi=300)
            plt.close()

            pcm = plt.pcolormesh(x_grid, t_light, np.real(psi), cmap=parula_map,vmin=-np.amax(np.real(psi)), vmax=np.amax(np.real(psi)), shading='auto')
            cbar = plt.colorbar(pcm, shrink=1)
            cbar.set_label('$\psi_R(x, t)$', rotation=0, size=20, labelpad=-27, y=1.1)
            plt.xlim([x_grid[0], x_grid[-1]])
            plt.xlabel('$x$', size='20')
            plt.ylabel('$t$', size='20')
            plt.grid(linestyle='--', alpha=0.5)
            plt.savefig(file + subfile + '/pdnlS_real_spacetime.png', dpi=300)
            plt.close()

            pcm = plt.pcolormesh(x_grid, t_light, np.imag(psi), cmap=parula_map, vmin=-np.amax(np.imag(psi)), vmax=np.amax(np.imag(psi)), shading='auto')
            cbar = plt.colorbar(pcm, shrink=1)
            cbar.set_label('$\psi_I(x, t)$', rotation=0, size=20, labelpad=-27, y=1.1)
            plt.xlim([x_grid[0], x_grid[-1]])
            plt.xlabel('$x$', size='20')
            plt.ylabel('$t$', size='20')
            plt.grid(linestyle='--', alpha=0.5)
            plt.savefig(file + subfile + '/pdnlS_imag_spacetime.png', dpi=300)
            plt.close()

            plt.plot(x_grid, np.real(psi)[-1, :], label="$\\textrm{Re}\ A$", color="b", zorder=5, lw=2)
            plt.plot(x_grid, np.imag(psi)[-1, :], label="$\\textrm{Im}\ A$", color="r", zorder=5, lw=2)
            plt.xlabel('$x$', size='25')
            plt.legend(fontsize=18)
            plt.xlim([-200, 200])
            plt.grid(alpha=0.3)
            plt.tight_layout()
            plt.savefig(file + subfile + '/final_profiles.png', dpi=200)
            plt.close()

            U_1_init = U1_light[-1, :]