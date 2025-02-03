from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':

    # Definiendo parámetros
    project_name = '/GPE'
    disc = 'D:/'
    route = 'mnustes_science/simulation_data/FD'
    eq = 'GPE'
    t_rate = 1


    # Definición de la grilla
    [tmin, tmax, dt] = [0, 5, 0.002]
    [xmin, xmax, dx] = [-20, 20, 0.05]
    t_grid = np.arange(tmin, tmax + dt, dt)
    x_grid = np.arange(xmin, xmax, dx)
    T = tmax
    Nt = t_grid.shape[0]
    Nx = x_grid.shape[0]

    alpha = 0.5
    beta = 1.0
    V_0 = 1.0
    mu = 1.0

    # Initial Conditions Pattern
    U_1_init = np.sqrt(2 * mu) / np.cosh(np.sqrt(2 * mu) * (x_grid - 2)) + 1j * 0.0
    #U_1_init = 1 / np.cosh(x_grid) + 1j * 0.0

    # Empaquetamiento de parametros, campos y derivadas para integración
    L = xmax - xmin
    D2 = sparse_DD_neumann(Nx, dx)
    operators = np.array([D2])
    fields_init = [U_1_init]
    grids = [t_grid, x_grid, 0]

    V = V_0 * (x_grid / 2) ** 2
    #plt.plot(x_grid, U_1_init)
    #plt.plot(x_grid, V)
    #plt.show()

    #delta_potential = np.exp(- (x_grid - 10) ** 2 / (2 * 6 ** 2)) - np.exp(- (x_grid + 10) ** 2 / (2 * 6 ** 2))
    parameters = [alpha, beta, V, mu]

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
    U1 = np.array(fields_history)[:, 0]
    U_complex = U1
    t_grid = time_grid
    modulo = np.absolute(U_complex)
    arg = np.angle(U_complex)
    arg = (2 * np.pi + arg) * (arg < 0) + arg * (arg > 0)

    lightness = 1
    U1_light = U1[0:-1:lightness]
    U_complex_light = U_complex[0:-1:lightness]
    t_light = time_grid[0:-1:lightness]
    modulo_light = modulo[0:-1:lightness]
    arg_light = arg[0:-1:lightness]

    # Guardando datos
    file = disc + route + project_name
    subfile = "/test"
    if not os.path.exists(file + subfile):
        os.makedirs(file + subfile)

    np.savetxt(file + subfile + '/field_real.txt', U1, delimiter=',')
    np.savetxt(file + subfile + '/X.txt', x_grid, delimiter=',')
    np.savetxt(file + subfile + '/T.txt', t_grid, delimiter=',')

    pcm = plt.pcolormesh(x_grid, t_light, modulo_light, cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$|A|$', rotation=0, size=20, labelpad=-27, y=1.1)
    plt.xlim([x_grid[0], x_grid[-1]])
    plt.xlabel('$x$', size='20')
    plt.ylabel('$t$', size='20')
    plt.grid(linestyle='--', alpha=0.5)
    plt.savefig(file + subfile + '/module_spacetime.png', dpi=300)
    plt.close()

    pcm = plt.pcolormesh(x_grid, t_light, arg_light, cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$\\textrm{arg}(A)$', rotation=0, size=20, labelpad=-20, y=1.1)
    plt.xlim([x_grid[0], x_grid[-1]])
    plt.xlabel('$x$', size='20')
    plt.ylabel('$t$', size='20')
    plt.grid(linestyle='--', alpha=0.5)
    plt.savefig(file + subfile + '/arg_spacetime.png', dpi=300)
    plt.close()