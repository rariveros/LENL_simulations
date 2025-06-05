from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':

    # Definiendo parámetros
    project_name = '/SH'
    disc = 'D:/'                                        # DISCO DE TRABAJO
    route = 'mnustes_science/simulation_data/FD'        # CARPETA DE TRABAJO
    eq = 'SH'                         # ECUACION
    t_rate = 10                                         # CADA CUANTAS ITERACIONES GUARDA
    dt = 0.002
    T = 600
    dx = 0.5 #en milimetros
    ies = [1]
    jotas = [1] #np.arange(0.25, 0.18, - 0.005)

    [tmin, tmax, dt] = [0, T, dt]
    [xmin, xmax, dx] = [-120, 120, dx]
    t_grid = np.arange(tmin, tmax + dt, dt)
    x_grid = np.arange(xmin, xmax, dx)
    T = tmax
    Nt = t_grid.shape[0]
    Nx = x_grid.shape[0]

    print("N° of simulations: " + str(len(ies) * len(jotas)))
    for i in ies:
        U_1_init = 0.01 * np.random.rand(Nx)
        for j in jotas:
            delta = 1.0
            gamma = 0.0
            q = 0.5
            epsilon_0 = 0.1
            dist = 50
            sigma = 16

            # Empaquetamiento de parametros, campos y derivadas para integración
            L = xmax - xmin
            D1 = sparse_D_neumann(Nx, dx)
            D2 = sparse_DD_neumann(Nx, dx)
            D4 = sparse_DDDD_neumann(Nx, dx)
            operators = np.array([D1, D2, D4])
            fields_init = [U_1_init]
            grids = [t_grid, x_grid, 0]
            phi = 0
            epsilon = epsilon_0 * (np.exp(- (x_grid - dist / 2) ** 2 / (2 * sigma ** 2)) + np.exp(1j * phi) * np.exp(- (x_grid + dist / 2) ** 2 / (2 * sigma ** 2)))

            parameters = [delta, gamma, q, epsilon]

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

            # Definiendo variables finales
            modulo_light = np.absolute(U1_light)

            # Guardando datos
            file = disc + route + project_name
            subfile = "/localized"
            parameters_np = np.array([delta, gamma, q, epsilon_0, dist, sigma])

            if not os.path.exists(file + subfile):
                os.makedirs(file + subfile)
            np.savetxt(file + subfile + '/field.txt', U1_light, delimiter=',')
            np.savetxt(file + subfile + '/parameters.txt', parameters_np, delimiter=',')
            np.savetxt(file + subfile + '/X.txt', x_grid, delimiter=',')
            np.savetxt(file + subfile + '/T.txt', t_light, delimiter=',')

            #################### Guardando Gráficos ####################
            pcm = plt.pcolormesh(x_grid, t_light, np.real(U1_light), cmap=parula_map, shading='auto')
            cbar = plt.colorbar(pcm, shrink=1)
            cbar.set_label('$|A|$', rotation=0, size=20, labelpad=-27, y=1.1)
            plt.xlim([x_grid[0], x_grid[-1]])
            plt.xlabel('$x$', size='20')
            plt.ylabel('$t$', size='20')
            plt.grid(linestyle='--', alpha=0.5)
            plt.savefig(file + subfile + '/A_module.png', dpi=300)
            plt.close()

            plt.plot(x_grid, np.real(U1_light[-1, :]), label="$u(x,t)$", color="b", zorder=5, lw=2)
            plt.xlabel('$x$', size='25')
            plt.legend(fontsize=18)
            plt.xlim([x_grid[0], x_grid[-1]])
            plt.grid(alpha=0.3)
            plt.tight_layout()
            plt.savefig(file + subfile + '/final_profiles.png', dpi=200)
            plt.close()

            U_1_init = U1_light[-1, :]