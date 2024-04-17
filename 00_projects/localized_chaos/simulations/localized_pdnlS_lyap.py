from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    #ti.init(arch=ti.gpu)
    # Definiendo parámetros
    project_name = '/localized_chaos'
    disc = 'C:/'
    route = 'mnustes_science/simulation_data/FD'
    eq = 'PDNLS'
    t_rate = 10
    [tmin, tmax, dt] = [0, 1000, 0.02]
    [xmin, xmax, dx] = [-60, 60, 0.2]
    t_grid = np.arange(tmin, tmax + dt, dt)
    x_grid = np.arange(xmin, xmax, dx)
    T = tmax
    Nt = t_grid.shape[0]
    Nx = x_grid.shape[0]
    D_kys = []
    param = []
    N_condit = int(1 * 2 * Nx)
    x_condit = np.arange(0, N_condit)/(2*Nx)
    U_1_init = 0.01 * np.random.rand(Nx)  # Z_r_0[-1, :]#
    U_2_init = 0.01 * np.random.rand(Nx)  # Z_i_0[-1, :]#Z_i_0[-1, :]#
    ies = [0.275, 0.300]#np.arange(0.375, 0.525, 0.025)#
    jotas = [8]
    m = 10
    i = 0
    for gamma_i in ies:
        for sigma_i in jotas:
            print("N° of simulations: " + str(len(ies)*len(jotas)))
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
            mu[0:10] = 10
            mu[-10:] = 10

            parameters = [alpha, beta, gamma, mu_0, nu]
            [alpha_str, beta_str, mu_str, nu_str, sigma_str, gamma_str] = pdnlS_str_parameters([alpha, beta, mu_0, nu, sigma, gamma_0], 0)

            # SIMULANDO Y CALCULANDO ESPECTRO
            now = datetime.datetime.now()
            print('Hora de Inicio: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
            time_init = time.time()
            print("Simulando y calculando espectro...")

            time_grid, lyapunov, fields_history, time_grid_fields = RK4_FD_lyapunov(eq, fields_init, parameters, grids, dt, Nt, operators, t_rate, N_condit, m)

            U1 = np.array(fields_history)[:, 0]
            U2 = np.array(fields_history)[:, 1]
            U_complex = U1 + 1j * U2
            modulo = np.absolute(U_complex)
            arg = np.angle(U_complex)
            arg = (2 * np.pi + arg) * (arg < 0) + arg * (arg > 0)

            lightness = 1
            U1_light = U1[::lightness]
            U2_light = U2[::lightness]
            U_complex_light = U_complex[::lightness]
            t_light = time_grid_fields[::lightness]
            modulo_light = modulo[::lightness]
            arg_light = arg[::lightness]

            average_Nt_i = int(0.3 * len(time_grid))
            lyapunov = np.array(lyapunov)
            lyap_mean = np.mean(lyapunov[average_Nt_i:], axis=0)

            sum = 0
            if lyap_mean[0] <= 0:
                D_ky = 0
            else:
                for i in range(len(x_condit)):
                    if sum >= 0:
                        sum = sum + lyap_mean[i]
                        first_neg = lyap_mean[i + 1]
                        p = i
                D_ky = p + sum / first_neg
            D_kys.append(D_ky)

            file = disc + route + project_name
            subfile = pdnlS_name([alpha, beta, mu_0, nu, sigma, gamma_0], "ABMNSG")
            parameters_np = np.array([alpha, beta, gamma_0, mu_0, nu, sigma])

            if m != 10:
                if not os.path.exists(file + subfile):
                    os.makedirs(file + subfile)
                print("Guardando datos...")
                np.savetxt(file + subfile + '/lyapunov_mean.txt', lyap_mean, delimiter=',')
                np.savetxt(file + subfile + '/lyapunov_espectrum.txt', lyapunov, delimiter=',')
                np.savetxt(file + subfile + '/D_ky.txt', np.array([gamma_0, sigma, D_ky]), delimiter=',')
                np.savetxt(file + subfile + '/field_real.txt', U1[-1, :], delimiter=',')
                np.savetxt(file + subfile + '/field_img.txt', U2[-1, :], delimiter=',')

                plt.plot(time_grid, lyapunov)
                plt.savefig(file + subfile + '/lyapunovs.png', dpi=200)
                plt.close()

                plt.plot(x_condit, np.flip(np.sort(lyap_mean)))
                plt.savefig(file + subfile + '/lyapunov_spectrum.png', dpi=200)
                plt.close()

                pcm = plt.pcolormesh(x_grid, t_light, modulo_light, cmap=parula_map, shading='auto')
                cbar = plt.colorbar(pcm, shrink=1)
                cbar.set_label('$|A|$', rotation=0, size=20, labelpad=-27, y=1.1)
                plt.xlim([x_grid[0], x_grid[-1]])
                plt.xlabel('$x$', size='20')
                plt.ylabel('$t$', size='20')
                plt.ylim([t_grid[-12000], t_grid[-1]])
                plt.grid(linestyle='--', alpha=0.5)
                plt.title(
                    '$\gamma = ' + str(gamma_str) + '$' + '    $\\nu = ' + str(nu_str) + '$' + '    $\mu = ' + str(mu_str) + '$',
                    size='20')
                plt.savefig(file + subfile + '/module_spacetime.png', dpi=100)
                plt.close()

                pcm = plt.pcolormesh(x_grid, time_grid_fields, modulo, cmap=parula_map, shading='auto')
                cbar = plt.colorbar(pcm, shrink=1)
                cbar.set_label('$|A|$', rotation=0, size=20, labelpad=-27, y=1.1)
                plt.xlim([x_grid[0], x_grid[-1]])
                plt.xlabel('$x$', size='20')
                plt.ylabel('$t$', size='20')
                plt.ylim([t_grid[-1000], t_grid[-1]])
                plt.grid(linestyle='--', alpha=0.5)
                plt.title(
                    '$\gamma = ' + str(gamma_str) + '$' + '    $\\nu = ' + str(nu_str) + '$' + '    $\mu = ' + str(mu_str) + '$',
                    size='20')
                plt.savefig(file + subfile + '/module_spacetime_short.png', dpi=100)
                plt.close()

            now = datetime.datetime.now()
            print('Hora de Término: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
            time_fin = time.time()
            print(str(time_fin - time_init) + ' seg')

            tmax = 1000
            t_grid = np.arange(tmin, tmax + dt, dt)
            Nt = t_grid.shape[0]

            U_1_init = U1_light[-1, :]
            U_2_init = U2_light[-1, :]
            m = 15000
    plt.scatter(param, D_kys)
    plt.savefig('D_ky_s12.png', dpi=200)
    plt.close()