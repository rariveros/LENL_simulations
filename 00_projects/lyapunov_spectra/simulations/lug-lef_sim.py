from functions import *
from back_process import *
from time_integrators import *
from directories_lyap import *

if __name__ == '__main__':

    # Definiendo parámetros
    project_name = main_directory + directoy_sim
    disc = 'E:/'
    route = 'mnustes_science/simulation_data/FD'
    eq = 'lug_lef'
    save_rate = 5
    [tmin, tmax, dt] = [0, 100, 0.005]
    [xmin, xmax, dx] = [0, 60, 0.1]

    powers = [7, 8, 9, 10]

    for power_i in powers:
        detuning = 1.1
        power = power_i
        C = 1

        # Definición de la grilla

        t_grid = np.arange(tmin, tmax + dt, dt)
        x_grid = np.arange(xmin, xmax, dx)
        T = tmax
        Nt = t_grid.shape[0]
        Nx = x_grid.shape[0]

        # Initial Conditions Pattern
        U_1_init = 0.1 * np.random.rand(Nx)
        U_2_init = 0.1 * np.random.rand(Nx)

        # Empaquetamiento de parametros, campos y derivadas para integración
        L = xmax - xmin
        D2 = sparse_DD(Nx, dx)
        operators = np.array([D2])
        fields_init = [U_1_init, U_2_init]
        grids = [t_grid, x_grid, 0]
        parameters = [detuning, power, C]

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

        # Guardando datos
        file = disc + route + project_name
        subfile = "/lug_lef/E=" + str(power_i)
        parameters_np = np.array([detuning, power])
        if not os.path.exists(file + subfile):
            os.makedirs(file + subfile)
        np.savetxt(file + subfile + '/field_real.txt', U1_light, delimiter=',')
        np.savetxt(file + subfile + '/field_img.txt', U2_light, delimiter=',')
        np.savetxt(file + subfile + '/parameters.txt', parameters_np, delimiter=',')
        np.savetxt(file + subfile + '/final_envelope.txt', amplitude_envelope_1, delimiter=',')
        np.savetxt(file + subfile + '/X.txt', x_grid, delimiter=',')
        np.savetxt(file + subfile + '/T.txt', t_light, delimiter=',')

        pcm = plt.pcolormesh(x_grid, t_light, modulo_light_1, cmap=parula_map, shading='auto')
        cbar = plt.colorbar(pcm, shrink=1)
        cbar.set_label('$|A|$', rotation=0, size=20, labelpad=-27, y=1.1)
        cbar.ax.tick_params(labelsize=15)
        plt.xlim([x_grid[0], x_grid[-1]])
        plt.xlabel('$x$', size='25')
        plt.ylabel('$t$', size='25')
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.grid(linestyle='--', alpha=0.2, color='k')
        plt.tight_layout()
        plt.savefig(file + subfile + '/module_spacetime.png', dpi=300)
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
        plt.savefig(file + subfile + '/arg_spacetime.png', dpi=300)
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
        plt.savefig(file + subfile + '/real_spacetime.png', dpi=300)
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
        plt.savefig(file + subfile + '/img_spacetime.png', dpi=300)
        plt.close()