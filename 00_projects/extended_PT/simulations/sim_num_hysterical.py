from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':

    # Definiendo parámetros
    project_name = '/PDNLS_extended_PT/extras_02'
    disc = 'C:/'                                        # DISCO DE TRABAJO
    route = 'mnustes_science/simulation_data/FD'        # CARPETA DE TRABAJO
    eq = 'PDNLS'                                        # ECUACION
    t_rate = 100                                        # CADA CUANTAS ITERACIONES GUARDA
    dt = 0.01
    T = 5000
    dx = 0.25 #en milimetros
    #directory = "C:/mnustes_science/simulation_data/FD/PDNLS_extended_PT/extras_02/dist_bifurcation/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.0180/sigma=6.000/gamma=0.1850/dist=30.2920"
    #Z_r_0 = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    #Z_i_0 = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    ies = np.arange(16, 80, 2)

    # Definición de la grilla
    [tmin, tmax, dt] = [0, T, dt]
    [xmin, xmax, dx] = [-100, 100, dx]
    t_grid = np.arange(tmin, tmax + dt, dt)
    x_grid = np.arange(xmin, xmax, dx)
    T = tmax
    Nt = t_grid.shape[0]
    Nx = x_grid.shape[0]
    U_1_init = 0.01 * np.random.rand(Nx)#
    U_2_init = 0.01 * np.random.rand(Nx)#
    print("N° of simulations: " + str(len(ies)))
    for dist in ies:
        alpha = 1  #5.721
        beta = 1
        nu = 0.32  #0.04812#0.0327449 #0.0052
        mu = 0.1
        gamma_0 = 0.28
        sigma = 3
        #dist = 20

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
        phi = np.pi
        gamma_l = 1
        gamma_r = 1
        gamma_real = gamma_0 * (gamma_l * np.exp(- (x_grid - dist / 2) ** 2 / (2 * sigma ** 2)) + gamma_r * np.cos(phi) * np.exp(- (x_grid + dist / 2) ** 2 / (2 * sigma ** 2)))
        gamma_img = gamma_r * np.sin(phi) * np.exp(- (x_grid + dist / 2) ** 2 / (2 * sigma ** 2))
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
        subfile = pdnlS_bigauss_name([alpha, beta, mu, nu, sigma, gamma_0, dist], "ABMNSGD")
        parameters_np = np.array([alpha, beta, gamma_0, dist,  mu, nu, sigma])

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
        plt.savefig(file + subfile + '/forcing_01.png', dpi=300)
        plt.close()

        pcm = plt.pcolormesh(x_grid, t_light, modulo_light_1, cmap=parula_map, shading='auto')
        cbar = plt.colorbar(pcm, shrink=1)
        cbar.set_label('$|A|$', rotation=0, size=20, labelpad=-27, y=1.1)
        plt.xlim([x_grid[0], x_grid[-1]])
        plt.xlabel('$x$', size='20')
        plt.ylabel('$t$', size='20')
        plt.grid(linestyle='--', alpha=0.5)
        plt.title('$\gamma_0 = ' + gamma_str + '\ \\alpha = ' + alpha_str  + '\ \\textrm{mm}^{2}' + '\ \\beta = ' + beta_str + '\ \\nu = ' + nu_str + '$', size='12')
        plt.savefig(file + subfile + '/module_spacetime.png', dpi=300)
        plt.close()

        pcm = plt.pcolormesh(x_grid, t_light, arg_light_1, cmap=parula_map, shading='auto')
        cbar = plt.colorbar(pcm, shrink=1)
        cbar.set_label('$\\textrm{arg}(A)$', rotation=0, size=20, labelpad=-20, y=1.1)
        plt.xlim([x_grid[0], x_grid[-1]])
        plt.xlabel('$x$', size='20')
        plt.ylabel('$t$', size='20')
        plt.grid(linestyle='--', alpha=0.5)
        plt.title('$\gamma_0 = ' + gamma_str + '\ \\alpha = ' + alpha_str  + '\ \\textrm{mm}^{2}' + '\ \\beta = ' + beta_str + '\ \\nu = ' + nu_str + '$', size='12')
        plt.savefig(file + subfile + '/arg_spacetime.png', dpi=300)
        plt.close()

        pcm = plt.pcolormesh(x_grid, t_light, U1_light, cmap=parula_map,vmin=-np.amax(U1_light), vmax=np.amax(U1_light), shading='auto')
        cbar = plt.colorbar(pcm, shrink=1)
        cbar.set_label('$A_R(x, t)$', rotation=0, size=20, labelpad=-27, y=1.1)
        plt.xlim([x_grid[0], x_grid[-1]])
        plt.xlabel('$x$', size='20')
        plt.ylabel('$t$', size='20')
        plt.grid(linestyle='--', alpha=0.5)
        plt.title('$\gamma_0 = ' + gamma_str + '\ \\alpha = ' + alpha_str  + '\ \\textrm{mm}^{2}' + '\ \\beta = ' + beta_str + '\ \\nu = ' + nu_str + '$', size='12')
        plt.savefig(file + subfile + '/real_spacetime.png', dpi=300)
        plt.close()

        pcm = plt.pcolormesh(x_grid, t_light, U2_light, cmap=parula_map, shading='auto')
        cbar = plt.colorbar(pcm, shrink=1)
        cbar.set_label('$A_I(x, t)$', rotation=0, size=20, labelpad=-27, y=1.1)
        plt.xlim([x_grid[0], x_grid[-1]])
        plt.xlabel('$x$', size='20')
        plt.ylabel('$t$', size='20')
        plt.grid(linestyle='--', alpha=0.5)
        plt.title('$\gamma_0 = ' + gamma_str + '\ \\alpha = ' + alpha_str  + '\ \\textrm{mm}^{2}' + '\ \\beta = ' + beta_str + '\ \\nu = ' + nu_str + '$', size='12')
        plt.savefig(file + subfile + '/img_spacetime.png', dpi=300)
        plt.close()
