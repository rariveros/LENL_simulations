from back_process import *
from functions import *
from time_integrators import *


if __name__ == '__main__':
    # Midiendo tiempo inicial
    now = datetime.datetime.now()
    print('Hora de Inicio: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
    time_init = time.time()

    # Definiendo parámetros
    project_name = '/C_nonresonant_article'
    disc = 'C:/'
    eq = 'PNDLS_envelope_nonresonant'
    for i in [0.30]:
        print('i = ' + str(i))
        alpha = 1
        beta = 1
        gamma = 0
        mu = 0.15
        nu = 0.4
        sigma = 10

        print('delta =' + str(gamma - (mu + (1 / sigma) * np.sqrt(nu / alpha))))

        gamma_str = str(int(gamma * 1000) * 0.001)
        nu_str = str(int(nu * 1000) * 0.001)
        mu_str = str(int(mu * 1000) * 0.001)

        print('gamma = ' + gamma_str)
        print('nu = ' + nu_str)
        print('mu = ' + mu_str)

        # Ploteo de Lengua de Arnold
        plot_parameters = 'no'
        if plot_parameters == 'si':
            arnold_tongue_show(gamma, mu, nu)

        # Definición de la grilla
        [tmin, tmax, dt] = [0, 1, 0.001]
        [xmin, xmax, dx] = [-60, 60, 0.1]
        t_grid = np.arange(tmin, tmax + dt, dt)
        x_grid = np.arange(xmin, xmax, dx)
        T = tmax
        Nt = t_grid.shape[0]
        Nx = x_grid.shape[0]

        # Initial Conditions
        U_1 = np.zeros((Nt, Nx))
        U_2 = np.zeros((Nt, Nx))

        U_10 = 0.01 * np.random.rand(Nx)
        U_20 = 0.01 * np.random.rand(Nx)
        U_1[0, :] = U_10
        U_2[0, :] = U_20

        # Empaquetamiento de parametros, campos y derivadas para integración
        fields = np.array([U_1, U_2])
        parameters = [alpha, beta, gamma, mu, nu, sigma]
        D1 = sparse_D(Nx, dx)
        D2 = sparse_DD(Nx, dx)
        operators = np.array([D1, D2])

        # Integración temporal
        final_fields = RK4_complexfields_FD(eq, fields, parameters, x_grid, dt, Nt, operators)

        # Midiendo tiempo final
        now = datetime.datetime.now()
        print('Hora de Término: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
        time_fin = time.time()
        print(str(time_fin - time_init) + ' seg')

        # Aligerando campos
        ratio = 500
        U1_light = final_fields[0, 0:-1:ratio, :]
        U2_light = final_fields[1, 0:-1:ratio, :]
        U_complex = U1_light + 1j * U2_light
        t_light = t_grid[0:-1:ratio]

        # Definiendo variables finales
        modulo_light = np.absolute(U_complex)
        arg_light = np.angle(U_complex)
        arg_light = (2*np.pi + arg_light) * (arg_light < 0) + arg_light*(arg_light > 0)
        analytical_signal = hilbert(U1_light[-1, :])
        amplitude_envelope = np.abs(analytical_signal)

        # Guardando datos
        file = disc + 'Users/mnustes_science/PT_fluids/mnustes_science/simulation_data/FD' + project_name
        subfile = nombre_pndls_gaussian(gamma, mu, nu, sigma)
        parameters_np = np.array([alpha, beta, gamma, mu, nu])
        if not os.path.exists(file + subfile):
            os.makedirs(file + subfile)
        np.savetxt(file + subfile + '/field_real.txt', U1_light, delimiter=',')
        np.savetxt(file + subfile + '/field_img.txt', U2_light, delimiter=',')
        np.savetxt(file + subfile + '/parameters.txt', parameters_np, delimiter=',')
        np.savetxt(file + subfile + '/final_envelope.txt', amplitude_envelope, delimiter=',')
        np.savetxt(file + subfile + '/X.txt', x_grid, delimiter=',')
        np.savetxt(file + subfile + '/T.txt', t_light, delimiter=',')

        # Gráficos
        plt.plot(x_grid, modulo_light[0, :], label='$R(x, 0)$')
        plt.plot(x_grid, arg_light[0, :], label='$\phi(x, 0)$')
        plt.plot(x_grid, U1_light[0, :], label='$\psi(x, 0)_{R}$')
        plt.plot(x_grid, U2_light[0, :], label='$\psi(x, 0)_{I}$')
        plt.legend()
        plt.grid()
        plt.xlim([x_grid[0], x_grid[-1]])
        #plt.savefig(file + subfile + '/initial_conditions.png', dpi=300)
        plt.close()

        pcm = plt.pcolormesh(x_grid, t_light, modulo_light, cmap='jet', shading='auto')
        cbar = plt.colorbar(pcm, shrink=1)
        cbar.set_label('$|C|$', rotation=0, size=20, labelpad=-27, y=1.1)
        plt.xlim([x_grid[0], x_grid[-1]])
        plt.xlabel('$x$', size='20')
        plt.ylabel('$t$', size='20')
        plt.grid(linestyle='--', alpha=0.5)
        plt.title('$\gamma = ' + str(gamma_str) + '$' + '    $\\nu = ' + str(nu_str) + '$' + '    $\mu = ' + str(mu_str) + '$', size='20')
        plt.savefig(file + subfile + '/module_spacetime.png', dpi=300)
        plt.close()

        pcm = plt.pcolormesh(x_grid, t_light, U1_light, cmap='seismic', shading='auto')
        cbar = plt.colorbar(pcm, shrink=1)
        cbar.set_label('$C_R(x, t)$', rotation=0, size=20, labelpad=-27, y=1.1)
        plt.xlim([x_grid[0], x_grid[-1]])
        plt.xlabel('$x$', size='20')
        plt.ylabel('$t$', size='20')
        plt.grid(linestyle='--', alpha=0.5)
        plt.title('$\gamma = ' + str(gamma_str) + '$' + '    $\\nu = ' + str(nu_str) + '$' + '    $\mu = ' + str(mu_str) + '$', size='20')
        plt.savefig(file + subfile + '/real_spacetime.png', dpi=300)
        plt.close()

        pcm = plt.pcolormesh(x_grid, t_light, arg_light, cmap='jet', shading='auto')
        cbar = plt.colorbar(pcm, shrink=1)
        cbar.set_label('$arg(C)$', rotation=0, size=20, labelpad=-20, y=1.1)
        plt.xlim([x_grid[0], x_grid[-1]])
        plt.xlabel('$x$', size='20')
        plt.ylabel('$t$', size='20')
        plt.grid(linestyle='--', alpha=0.5)
        plt.title('$\gamma = ' + str(gamma_str) + '$' + '    $\\nu = ' + str(nu_str) + '$' + '    $\mu = ' + str(mu_str) + '$', size='20')
        plt.savefig(file + subfile + '/arg_spacetime.png', dpi=300)
        plt.close()

        nu_positive_grid = np.arange(0, 2, 0.01)
        nu_negative_grid = - np.flip(nu_positive_grid)
        nu_grid = np.append(nu_negative_grid, nu_positive_grid)
        plt.plot(nu_positive_grid, np.sqrt(nu_positive_grid ** 2 + mu ** 2), c='k', linestyle='--')
        plt.fill_between(nu_positive_grid, np.ones(len(nu_positive_grid)) * mu, np.sqrt(nu_positive_grid ** 2 + mu ** 2),
                         facecolor=(92 / 255, 43 / 255, 228 / 255, 0.4))
        plt.plot(nu_negative_grid, np.sqrt(nu_negative_grid ** 2 + mu ** 2), c='k', linestyle='--')
        plt.fill_between(nu_negative_grid, np.ones(len(nu_negative_grid)) * mu, np.sqrt(nu_negative_grid ** 2 + mu ** 2),
                         facecolor=(0, 1, 0, 0.4))
        plt.plot(nu_grid, np.ones(len(nu_grid)) * mu, c='k', linestyle='--')
        plt.fill_between(nu_grid, 2, np.sqrt(nu_grid ** 2 + mu ** 2),
                         facecolor=(1, 0, 0, 0.4))
        plt.fill_between(nu_grid, np.ones(len(nu_grid)) * mu, 0,
                         facecolor=(1, 1, 0, 0.4))
        plt.scatter(nu, gamma, c='k', zorder=10)
        plt.title('Arnold Tongue', size='25')
        plt.xlabel('$\\nu$', size='25')
        plt.ylabel('$\gamma$', size='25')
        plt.xlim([-1, 1])
        plt.ylim([0, 1])
        plt.grid(linestyle='--', alpha=0.5)
        plt.savefig(file + subfile + '/arnold_tongue.png')
        plt.close()