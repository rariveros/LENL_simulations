from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':

    # Definiendo parámetros
    project_name = '/PDNLS_bigaussian'
    disc = 'D:/'                                        # DISCO DE TRABAJO
    route = 'mnustes_science/simulation_data/FD'        # CARPETA DE TRABAJO
    eq = 'PDNLS'                                        # ECUACION
    t_rate = 100                                        # CADA CUANTAS ITERACIONES GUARDA
    dt = 0.01
    T = 1000
    dx = 1 #en milimetros

    pist_dist = [1, 2]                 # PISTONES DE DISTANCIA ENTRE GAUSSIANAS
    f_i = np.arange(15.18, 15.2, 0.001)                                      # FRECUENCIA DE INYECCION
    a = [14.5]                                             # AMPLITUD ANGULAR
    for i in range(len(f_i)):
        for j in range(len(a)):
            for pist_dist_i in pist_dist:
                alpha, beta, nu, gamma_0 = fluid_pdnls_parameters(f_i[i], a[j], d=20)         # PARAMETROS DE PDNLS DESDE FLUIDO
                beta = 1                                                                # RE-DEFINIR BETA COMO UNO
                mu = 0.10                                                               # PARAMETRO FENOMENOLOGICO
                piston_spac = 18.46                                                     # DISTANCIA ENTRE PISTONES
                n_pist = 1                                                              # CANTIDAD DE PISTONES ENCENDIDOS
                FWHM = n_pist * piston_spac
                sigma = FWHM / (2 * np.sqrt(2 * np.log(2)))
                print(sigma)# SE ASUME QUE LOS PISTONES USADOS CORRESPONDEN AL HALF WIDTH HALF MAXIMUM (HWHM)
                dist = pist_dist_i * piston_spac                # DISTANCIA ENTRE CENTROS DE GAUSSIANAS

                alpha_str = str(alpha).split(".")[0] + '.' + str(alpha).split(".")[1][0:3]
                beta_str = str(beta)
                nu_str = str(nu).split(".")[0] + '.' + str(nu).split(".")[1][0:3]
                gamma_str = str(gamma_0).split(".")[0] + '.' + str(gamma_0).split(".")[1][0:3]
                sigma_str = str(sigma).split(".")[0] + '.' + str(sigma).split(".")[1][0:1]

                # Definición de la grilla
                [tmin, tmax, dt] = [0, T, dt]
                [xmin, xmax, dx] = [-200, 200, dx]
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
                D2 = sparse_DD_neumann(Nx, dx)
                operators = np.array([D2])
                fields_init = [U_1_init, U_2_init]
                grids = [t_grid, x_grid, 0]
                gamma_real = gamma_0 * (np.exp(- (x_grid - dist / 2) ** 2 / (2 * sigma ** 2)) - np.exp(- (x_grid + dist / 2) ** 2 / (2 * sigma ** 2)))
                gamma_img = 0
                gamma = [gamma_real, 0]

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
                subfile = pdnlS_bigauss_name([n_pist, f_i[i], a[j], pist_dist_i], "exp_sim")
                parameters_np = np.array([alpha, beta, gamma_0, mu, nu, sigma])
                if not os.path.exists(file + subfile):
                    os.makedirs(file + subfile)
                np.savetxt(file + subfile + '/field_real.txt', U1_light, delimiter=',')
                np.savetxt(file + subfile + '/field_img.txt', U2_light, delimiter=',')
                np.savetxt(file + subfile + '/parameters.txt', parameters_np, delimiter=',')
                np.savetxt(file + subfile + '/final_envelope.txt', amplitude_envelope_1, delimiter=',')
                np.savetxt(file + subfile + '/X.txt', x_grid, delimiter=',')
                np.savetxt(file + subfile + '/T.txt', t_light, delimiter=',')

                # Guardando gráficos
                plt.plot(x_grid, gamma_real, color="k")
                plt.xlabel('$x\ \\textrm{(mm)}$', size='25')
                plt.xlim([x_grid[0], x_grid[-1]])
                plt.ylabel('$\gamma(x)$', size='25')
                plt.ylim([-gamma_0 * 1.1, gamma_0 * 1.1])
                plt.grid(linestyle='--', alpha=0.5)
                plt.savefig(file + subfile + '/forcing.png', dpi=300)
                plt.close()

                pcm = plt.pcolormesh(x_grid, t_light, modulo_light_1, cmap=parula_map, shading='auto')
                cbar = plt.colorbar(pcm, shrink=1)
                cbar.set_label('$|A|$', rotation=0, size=20, labelpad=-27, y=1.1)
                plt.xlim([x_grid[0], x_grid[-1]])
                plt.xlabel('$x$', size='20')
                plt.ylabel('$t$', size='20')
                plt.grid(linestyle='--', alpha=0.5)
                plt.title('$\gamma_0 = ' + gamma_str + '\ f_i = ' + str(f_i[i]) + '\ \\textrm{Hz}' +
                          '\ \\alpha = ' + alpha_str + '\ \\textrm{mm}^{2}' + '\ \\beta = ' + beta_str + '\ \\nu = ' + nu_str + '$',
                          size='10')
                plt.savefig(file + subfile + '/module_spacetime.png', dpi=300)
                plt.close()

                pcm = plt.pcolormesh(x_grid, t_light, arg_light_1, cmap=parula_map, shading='auto')
                cbar = plt.colorbar(pcm, shrink=1)
                cbar.set_label('$\\textrm{arg}(A)$', rotation=0, size=20, labelpad=-20, y=1.1)
                plt.xlim([x_grid[0], x_grid[-1]])
                plt.xlabel('$x$', size='20')
                plt.ylabel('$t$', size='20')
                plt.grid(linestyle='--', alpha=0.5)
                plt.title('$\gamma_0 = ' + gamma_str + '\ f_i = ' + str(f_i[i]) + '\ \\textrm{Hz}' +
                          '\ \\alpha = ' + alpha_str + '\ \\textrm{mm}^{2}' + '\ \\beta = ' + beta_str + '\ \\nu = ' + nu_str + '$',
                          size='10')
                plt.savefig(file + subfile + '/arg_spacetime.png', dpi=300)
                plt.close()

                pcm = plt.pcolormesh(x_grid, t_light, U1_light, cmap=parula_map,vmin=-np.amax(U1_light), vmax=np.amax(U1_light), shading='auto')
                cbar = plt.colorbar(pcm, shrink=1)
                cbar.set_label('$A_R(x, t)$', rotation=0, size=20, labelpad=-27, y=1.1)
                plt.xlim([x_grid[0], x_grid[-1]])
                plt.xlabel('$x$', size='20')
                plt.ylabel('$t$', size='20')
                plt.grid(linestyle='--', alpha=0.5)
                plt.title('$\gamma_0 = ' + gamma_str + '\ f_i = ' + str(f_i[i]) + '\ \\textrm{Hz}' +
                          '\ \\alpha = ' + alpha_str + '\ \\textrm{mm}^{2}' + '\ \\beta = ' + beta_str + '\ \\nu = ' + nu_str + '$',
                          size='10')
                plt.savefig(file + subfile + '/real_spacetime.png', dpi=300)
                plt.close()

                pcm = plt.pcolormesh(x_grid, t_light, U2_light, cmap=parula_map, shading='auto')
                cbar = plt.colorbar(pcm, shrink=1)
                cbar.set_label('$A_I(x, t)$', rotation=0, size=20, labelpad=-27, y=1.1)
                plt.xlim([x_grid[0], x_grid[-1]])
                plt.xlabel('$x$', size='20')
                plt.ylabel('$t$', size='20')
                plt.grid(linestyle='--', alpha=0.5)
                plt.title('$\gamma_0 = ' + gamma_str + '\ f_i = ' + str(f_i[i]) + '\ \\textrm{Hz}' +
                          '\ \\alpha = ' + alpha_str + '\ \\textrm{mm}^{2}' + '\ \\beta = ' + beta_str + '\ \\nu = ' + nu_str + '$',
                          size='10')
                plt.savefig(file + subfile + '/img_spacetime.png', dpi=300)
                plt.close()