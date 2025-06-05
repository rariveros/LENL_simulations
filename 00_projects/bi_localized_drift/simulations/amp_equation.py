import numpy as np

from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':

    # Definiendo parámetros
    project_name = '/bi_localized_drift/amplitude_eq/05_24_badone'
    disc = 'D:/'                                        # DISCO DE TRABAJO
    route = 'mnustes_science/simulation_data/FD'        # CARPETA DE TRABAJO
    eq = 'amplitude_bigaussian'                         # ECUACION
    t_rate = 200                                         # CADA CUANTAS ITERACIONES GUARDA
    dt = 0.02
    T = 5000
    dx = 1 #en milimetros
    ies = [80]#np.arange(80, 100, 2)
    jotas = [1] #np.arange(0.25, 0.18, - 0.005)

    [tmin, tmax, dt] = [0, T, dt]
    [xmin, xmax, dx] = [-130, 130, dx]
    t_grid = np.arange(tmin, tmax + dt, dt)
    x_grid = np.arange(xmin, xmax, dx)
    T = tmax
    Nt = t_grid.shape[0]
    Nx = x_grid.shape[0]
    #IC_directory = "D:/mnustes_science/simulation_data/FD/bi_localized_drift/amplitude_eq/05_24_dist/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.1000/sigma=16.000/gamma=0.2000/dist=78.000"
    #IC_R = np.loadtxt(IC_directory + "/ampd_field_real.txt", delimiter=',', dtype=complex)
    #IC_I = np.loadtxt(IC_directory + '/ampd_field_img.txt', delimiter=',', dtype=complex)

    #U_1_init = IC_R[-1, :]
    #U_2_init = IC_I[-1, :]
    print("N° of simulations: " + str(len(ies) * len(jotas)))
    U_1_init = 0.01 * (np.random.rand(Nx) + 1j * np.random.rand(Nx))
    U_2_init = 0.01 * (np.random.rand(Nx) + 1j * np.random.rand(Nx))
    for i in ies:
        for j in jotas:
            alpha = 6.524  #5.721
            beta = 1
            nu = 0.1 #0.014 #0.04812#0.0327449 #0.0052
            mu = 0.1
            gamma_0 = 0.24
            dist = i
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
            fields_init = [U_1_init, U_2_init]
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
            U2_light = np.array(fields_history)[:, 1]
            t_light = time_grid

            DU1 = []
            DU2 = []
            for i in range(len(t_light)):
                DU1.append(np.append(np.diff(U1_light[i, :]), 0) / dx)
                DU2.append(np.append(np.diff(U2_light[i, :]), 0) / dx)
            DU1 = np.array(DU1)
            DU2 = np.array(DU2)

            K = np.sqrt(nu / alpha)
            psiAR_1 = U1_light * np.exp(1j * K * (x_grid + dist / 2)) + np.conjugate(U1_light * np.exp(1j * K * (x_grid + dist / 2)))
            psiAR_3 = (1 / (8 * nu)) * U1_light ** 3 * np.exp(3 * 1j * K * (x_grid + dist / 2)) + np.conjugate((1 / (8 * nu)) * U1_light ** 3 * np.exp(3 * 1j * K * (x_grid + dist / 2)))
            psiAI_1 = U1_light * 0
            psiAI_3 = -(1 / (2 * mu)) * (2 * 1j * K * DU1 + 3 * np.abs(U1_light) ** 2 * U1_light) * np.exp(1j * K * (x_grid + dist / 2)) + np.conjugate(-(1 / (2 * mu)) * (2 * 1j * K * DU1 + 3 * np.abs(U1_light) ** 2 * U1_light) * np.exp(1j * K * (x_grid + dist / 2)))

            psiBI_1 = U2_light * np.exp(1j * K * (x_grid + dist / 2)) + np.conjugate(U2_light * np.exp(1j * K * (x_grid + dist / 2)))
            psiBI_3 = -(1 / (8 * nu)) * U2_light ** 3 * np.exp(3 * 1j * K * (x_grid + dist / 2)) + np.conjugate(- (1 / (8 * nu)) * U2_light ** 3 * np.exp(3 * 1j * K * (x_grid + dist / 2)))
            psiBR_1 = U2_light * 0
            psiBR_3 = -(1 / (2 * mu)) * (-2 * 1j * K * DU2 - 3 * np.abs(U2_light) ** 2 * U2_light) * np.exp(1j * K * (x_grid + dist / 2)) + np.conjugate(-(1 / (2 * mu)) * (-2 * 1j * K * DU2 - 3 * np.abs(U2_light) ** 2 * U2_light) * np.exp(1j * K * (x_grid + dist / 2)))

            psi_A = (psiAR_1 + psiAR_3) + 1j * (psiAI_1 + psiAI_3)
            psi_B = (psiBR_1 + psiBR_3) + 1j * (psiBI_1 + psiBI_3)
            psi = psi_A + psi_B

            # Definiendo variables finales
            modulo_light_1 = np.absolute(U1_light)
            modulo_light_2 = np.absolute(U2_light)
            arg_light_1 = np.angle(U1_light)
            arg_light_1 = (2 * np.pi + arg_light_1) * (arg_light_1 < 0) + arg_light_1 * (arg_light_1 > 0)
            arg_light_2 = np.angle(U2_light)
            arg_light_2 = (2 * np.pi + arg_light_2) * (arg_light_2 < 0) + arg_light_2 * (arg_light_1 > 0)

            # Guardando datos
            file = disc + route + project_name
            subfile = pdnlS_bigauss_name([alpha, beta, mu, nu, sigma, gamma_0, dist], "ABMNSGD")
            parameters_np = np.array([alpha, beta, gamma_0, dist, sigma,  mu, nu])

            if not os.path.exists(file + subfile):
                os.makedirs(file + subfile)
            np.savetxt(file + subfile + '/ampd_field_real.txt', U1_light, delimiter=',')
            np.savetxt(file + subfile + '/ampd_field_img.txt', U2_light, delimiter=',')
            np.savetxt(file + subfile + '/field_real.txt', np.real(psi), delimiter=',')
            np.savetxt(file + subfile + '/field_img.txt', np.imag(psi), delimiter=',')
            np.savetxt(file + subfile + '/parameters.txt', parameters_np, delimiter=',')
            np.savetxt(file + subfile + '/X.txt', x_grid, delimiter=',')
            np.savetxt(file + subfile + '/T.txt', t_light, delimiter=',')

            #################### Guardando Gráficos ####################
            pcm = plt.pcolormesh(x_grid, t_light, modulo_light_1 + modulo_light_2, cmap=parula_map, shading='auto')
            cbar = plt.colorbar(pcm, shrink=1)
            cbar.set_label('$|A|$', rotation=0, size=20, labelpad=-27, y=1.1)
            plt.xlim([x_grid[0], x_grid[-1]])
            plt.xlabel('$x$', size='20')
            plt.ylabel('$t$', size='20')
            plt.grid(linestyle='--', alpha=0.5)
            plt.savefig(file + subfile + '/A_module.png', dpi=300)
            plt.close()

            """
            pcm = plt.pcolormesh(x_grid, t_light, arg_light_1, cmap=parula_map, shading='auto')
            cbar = plt.colorbar(pcm, shrink=1)
            cbar.set_label('$\\textrm{arg}(A)$', rotation=0, size=20, labelpad=-20, y=1.1)
            plt.xlim([x_grid[0], x_grid[-1]])
            plt.xlabel('$x$', size='20')
            plt.ylabel('$t$', size='20')
            plt.grid(linestyle='--', alpha=0.5)
            plt.savefig(file + subfile + '/A_arg.png', dpi=300)
            plt.close()

            pcm = plt.pcolormesh(x_grid, t_light, np.real(U1_light), cmap=parula_map,vmin=-np.amax(np.real(U1_light)), vmax=np.amax(np.real(U1_light)), shading='auto')
            cbar = plt.colorbar(pcm, shrink=1)
            cbar.set_label('$A_R(x, t)$', rotation=0, size=20, labelpad=-27, y=1.1)
            plt.xlim([x_grid[0], x_grid[-1]])
            plt.xlabel('$x$', size='20')
            plt.ylabel('$t$', size='20')
            plt.grid(linestyle='--', alpha=0.5)
            plt.savefig(file + subfile + '/A_real.png', dpi=300)
            plt.close()

            pcm = plt.pcolormesh(x_grid, t_light, np.imag(U1_light), cmap=parula_map, shading='auto')
            cbar = plt.colorbar(pcm, shrink=1)
            cbar.set_label('$A_I(x, t)$', rotation=0, size=20, labelpad=-27, y=1.1)
            plt.xlim([x_grid[0], x_grid[-1]])
            plt.xlabel('$x$', size='20')
            plt.ylabel('$t$', size='20')
            plt.grid(linestyle='--', alpha=0.5)
            plt.savefig(file + subfile + '/A_img.png', dpi=300)
            plt.close()

            #####################################################

            pcm = plt.pcolormesh(x_grid, t_light, modulo_light_2, cmap=parula_map, shading='auto')
            cbar = plt.colorbar(pcm, shrink=1)
            cbar.set_label('$|A|$', rotation=0, size=20, labelpad=-27, y=1.1)
            plt.xlim([x_grid[0], x_grid[-1]])
            plt.xlabel('$x$', size='20')
            plt.ylabel('$t$', size='20')
            plt.grid(linestyle='--', alpha=0.5)
            plt.savefig(file + subfile + '/B_module.png', dpi=300)
            plt.close()

            pcm = plt.pcolormesh(x_grid, t_light, arg_light_2, cmap=parula_map, shading='auto')
            cbar = plt.colorbar(pcm, shrink=1)
            cbar.set_label('$\\textrm{arg}(A)$', rotation=0, size=20, labelpad=-20, y=1.1)
            plt.xlim([x_grid[0], x_grid[-1]])
            plt.xlabel('$x$', size='20')
            plt.ylabel('$t$', size='20')
            plt.grid(linestyle='--', alpha=0.5)
            plt.savefig(file + subfile + '/B_arg.png', dpi=300)
            plt.close()

            pcm = plt.pcolormesh(x_grid, t_light, np.real(U2_light), cmap=parula_map,vmin=-np.amax(np.real(U1_light)), vmax=np.amax(np.real(U1_light)), shading='auto')
            cbar = plt.colorbar(pcm, shrink=1)
            cbar.set_label('$A_R(x, t)$', rotation=0, size=20, labelpad=-27, y=1.1)
            plt.xlim([x_grid[0], x_grid[-1]])
            plt.xlabel('$x$', size='20')
            plt.ylabel('$t$', size='20')
            plt.grid(linestyle='--', alpha=0.5)
            plt.savefig(file + subfile + '/B_real.png', dpi=300)
            plt.close()

            pcm = plt.pcolormesh(x_grid, t_light, np.imag(U2_light), cmap=parula_map, shading='auto')
            cbar = plt.colorbar(pcm, shrink=1)
            cbar.set_label('$A_I(x, t)$', rotation=0, size=20, labelpad=-27, y=1.1)
            plt.xlim([x_grid[0], x_grid[-1]])
            plt.xlabel('$x$', size='20')
            plt.ylabel('$t$', size='20')
            plt.grid(linestyle='--', alpha=0.5)
            plt.savefig(file + subfile + '/B_img.png', dpi=300)
            plt.close()
            """
            #####################################################

            pcm = plt.pcolormesh(x_grid, t_light, np.abs(psi), cmap=parula_map, shading='auto')
            cbar = plt.colorbar(pcm, shrink=1)
            cbar.set_label('$|\psi|$', rotation=0, size=20, labelpad=-27, y=1.1)
            plt.xlim([x_grid[0], x_grid[-1]])
            plt.xlabel('$x$', size='20')
            plt.ylabel('$t$', size='20')
            plt.grid(linestyle='--', alpha=0.5)
            plt.savefig(file + subfile + '/pdnlS_module_spacetime.png', dpi=300)
            plt.close()

            """
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
            """
            plt.plot(x_grid, np.real(psi)[-1, :], label="$\\textrm{Re}\ A$", color="b", zorder=5, lw=2)
            plt.plot(x_grid, np.imag(psi)[-1, :], label="$\\textrm{Im}\ A$", color="r", zorder=5, lw=2)
            plt.xlabel('$x$', size='25')
            plt.legend(fontsize=18)
            plt.xlim([x_grid[0], x_grid[-1]])
            plt.grid(alpha=0.3)
            plt.tight_layout()
            plt.savefig(file + subfile + '/final_profiles.png', dpi=200)
            plt.close()

            U_1_init = U1_light[-1, :]
            U_2_init = U2_light[-1, :]