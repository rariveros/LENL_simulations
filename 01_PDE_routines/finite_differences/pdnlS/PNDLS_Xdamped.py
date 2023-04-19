from back_process import *
from functions import *
from time_integrators import *


from functions import *
from back_process import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':

    # Definiendo parámetros
    project_name = '/PNDLS_Xdamp'
    disc = 'C:/'
    eq = 'PNDLS_gaussian_Xdamp'

    gammas = [0.2]
    for i in gammas:
        alpha = 1
        beta = 1
        gamma_0 = i
        mu_01 = 0.1
        mu_02 = 0.01
        nu = 0.2
        sigma = 10

        gamma_str = str(int(gamma_0 * 1000) * 0.001)
        nu_str = str(int(nu * 1000) * 0.001)
        mu_str = str(int(mu_01 * 1000) * 0.001)
        print('gamma = ' + gamma_str)
        print('nu = ' + nu_str)
        print('mu = ' + mu_str)

        # Definición de la grilla
        [tmin, tmax, dt] = [0, 200, 0.001]
        [xmin, xmax, dx] = [-60, 60, 0.2]
        t_grid = np.arange(tmin, tmax + dt, dt)
        x_grid = np.arange(xmin, xmax, dx)
        T = tmax
        Nt = t_grid.shape[0]
        Nx = x_grid.shape[0]

        # Initial Conditions
        U_1_init = U_2_init = 0.01 * np.random.rand(Nx)

        # Empaquetamiento de parametros, campos y derivadas para integración
        L = xmax - xmin
        D2 = sparse_DD_neumann(Nx, dx)
        operators = np.array([D2])
        fields_init = [U_1_init, U_2_init]
        grids = [t_grid, x_grid, 0]
        A = i
        A_str = str(int(A * 1000) * 0.001)
        print('A = ' + A_str)
        F_sym = gamma_0 * np.exp(- x_grid ** 2 / (2 * sigma ** 2))
        F_1 = np.exp(-(x_grid - mu_01) ** 2 / (2 * sigma ** 2))
        F_2 = special.erfc(A * (x_grid - mu_01) / (np.sqrt(2) * sigma))
        F_asym = F_1 * F_2
        F_asym = gamma_0 * (1 / np.amax(F_asym)) * F_asym
        gamma = [F_sym, 0]

        parameters = [alpha, beta, gamma, mu_01, mu_02, nu]

        # Midiendo tiempo inicial
        now = datetime.datetime.now()
        print('Hora de Inicio: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
        time_init = time.time()

        final_fields, fields_history, time_grid = RK4_FD(eq, fields_init, parameters, grids, dt, Nt, operators)

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
        file = disc + 'Users/mnustes_science/PT_fluids/mnustes_science/simulation_data/FD' + project_name
        subfile = nombre_pndls_gaussian(gamma_0, mu_01, nu, sigma) #+ '/a=' + A_str
        parameters_np = np.array([alpha, beta, gamma_0, mu_01, nu])
        if not os.path.exists(file + subfile):
            os.makedirs(file + subfile)
        np.savetxt(file + subfile + '/field_real.txt', U1_light, delimiter=',')
        np.savetxt(file + subfile + '/field_img.txt', U2_light, delimiter=',')
        np.savetxt(file + subfile + '/parameters.txt', parameters_np, delimiter=',')
        np.savetxt(file + subfile + '/final_envelope.txt', amplitude_envelope_1, delimiter=',')
        np.savetxt(file + subfile + '/X.txt', x_grid, delimiter=',')
        np.savetxt(file + subfile + '/T.txt', t_light, delimiter=',')

        nu_positive_grid = np.arange(0, 2, 0.01)
        nu_negative_grid = - np.flip(nu_positive_grid)
        nu_grid = np.append(nu_negative_grid, nu_positive_grid)
        plt.plot(nu_positive_grid, np.sqrt(nu_positive_grid ** 2 + mu_01 ** 2), c='k', linestyle='--')
        plt.fill_between(nu_positive_grid, np.ones(len(nu_positive_grid)) * mu_01,
                         np.sqrt(nu_positive_grid ** 2 + mu_01 ** 2),
                         facecolor=(92 / 255, 43 / 255, 228 / 255, 0.4))
        plt.plot(nu_negative_grid, np.sqrt(nu_negative_grid ** 2 + mu_01 ** 2), c='k', linestyle='--')
        plt.fill_between(nu_negative_grid, np.ones(len(nu_negative_grid)) * mu_01,
                         np.sqrt(nu_negative_grid ** 2 + mu_01 ** 2),
                         facecolor=(0, 1, 0, 0.4))
        plt.plot(nu_grid, np.ones(len(nu_grid)) * mu_01, c='k', linestyle='--')
        plt.fill_between(nu_grid, 2, np.sqrt(nu_grid ** 2 + mu_01 ** 2),
                         facecolor=(1, 0, 0, 0.4))
        plt.fill_between(nu_grid, np.ones(len(nu_grid)) * mu_01, 0,
                         facecolor=(1, 1, 0, 0.4))
        plt.scatter(nu, gamma_0, c='k', zorder=10)
        plt.title('Arnold Tongue', size='25')
        plt.xlabel('$\\nu$', size='25')
        plt.ylabel('$\gamma$', size='25')
        plt.xlim([-1, 1])
        plt.ylim([0, 1])
        plt.grid(linestyle='--', alpha=0.5)
        plt.savefig(file + subfile + '/arnold_tongue.png')
        plt.close()

        pcm = plt.pcolormesh(x_grid, t_light, modulo_light_1, cmap=parula_map, shading='auto')
        cbar = plt.colorbar(pcm, shrink=1)
        cbar.set_label('$|A|$', rotation=0, size=20, labelpad=-27, y=1.1)
        plt.xlim([x_grid[0], x_grid[-1]])
        plt.xlabel('$x$', size='20')
        plt.ylabel('$t$', size='20')
        plt.grid(linestyle='--', alpha=0.5)
        plt.title(
            '$\gamma = ' + str(gamma_str) + '$' + '    $\\nu = ' + str(nu_str) + '$' + '    $\mu = ' + str(mu_str) + '$',
            size='20')
        plt.savefig(file + subfile + '/module_spacetime_01.png', dpi=300)
        plt.close()

        pcm = plt.pcolormesh(x_grid, t_light, arg_light_1, cmap=parula_map, shading='auto')
        cbar = plt.colorbar(pcm, shrink=1)
        cbar.set_label('$\\textrm{arg}(A)$', rotation=0, size=20, labelpad=-20, y=1.1)
        plt.xlim([x_grid[0], x_grid[-1]])
        plt.xlabel('$x$', size='20')
        plt.ylabel('$t$', size='20')
        plt.grid(linestyle='--', alpha=0.5)
        plt.title(
            '$\gamma = ' + str(gamma_str) + '$' + '    $\\nu = ' + str(nu_str) + '$' + '    $\mu = ' + str(mu_str) + '$',
            size='20')
        plt.savefig(file + subfile + '/arg_spacetime_01.png', dpi=300)
        plt.close()

        pcm = plt.pcolormesh(x_grid, t_light, U1_light, cmap=parula_map, shading='auto')
        cbar = plt.colorbar(pcm, shrink=1)
        cbar.set_label('$A_R(x, t)$', rotation=0, size=20, labelpad=-27, y=1.1)
        plt.xlim([x_grid[0], x_grid[-1]])
        plt.xlabel('$x$', size='20')
        plt.ylabel('$t$', size='20')
        plt.grid(linestyle='--', alpha=0.5)
        plt.title(
            '$\gamma = ' + str(gamma_str) + '$' + '    $\\nu = ' + str(nu_str) + '$' + '    $\mu = ' + str(mu_str) + '$',
            size='20')
        plt.savefig(file + subfile + '/real_spacetime_01.png', dpi=300)
        plt.close()

        pcm = plt.pcolormesh(x_grid, t_light, U2_light, cmap=parula_map, shading='auto')
        cbar = plt.colorbar(pcm, shrink=1)
        cbar.set_label('$A_I(x, t)$', rotation=0, size=20, labelpad=-27, y=1.1)
        plt.xlim([x_grid[0], x_grid[-1]])
        plt.xlabel('$x$', size='20')
        plt.ylabel('$t$', size='20')
        plt.grid(linestyle='--', alpha=0.5)
        plt.title(
            '$\gamma = ' + str(gamma_str) + '$' + '    $\\nu = ' + str(nu_str) + '$' + '    $\mu = ' + str(mu_str) + '$',
            size='20')
        plt.savefig(file + subfile + '/img_spacetime_01.png', dpi=300)
        plt.close()