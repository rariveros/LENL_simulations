from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':

    # Definiendo parámetros
    project_name = '/KdV'
    disc = 'C:/'
    eq = 'kdV'

    mu = 1

    # Definición de la grilla
    [tmin, tmax, dt] = [0, 30, 0.001]
    [xmin, xmax, dx] = [-30, 30, 0.1]
    t_grid = np.arange(tmin, tmax + dt, dt)
    x_grid = np.arange(xmin, xmax, dx)
    T = tmax
    Nt = t_grid.shape[0]
    Nx = x_grid.shape[0]

    # Initial Conditions
    a_1 = 1
    b_1 = -10
    U_1 = a_1 / (2 * np.cosh(0.5 * a_1 ** 0.5 * (x_grid - b_1)) ** 2)

    a_2 = 2
    b_2 = -20
    U_2 = a_2 / (2 * np.cosh(0.5 * a_2 ** 0.5 * (x_grid - b_2)) ** 2)

    U_init = U_1 + U_2
    # Empaquetamiento de parametros, campos y derivadas para integración
    L = xmax - xmin
    D1 = sparse_D(Nx, dx)
    D3 = sparse_DDD(Nx, dx)
    operators = np.array([D1, D3])
    fields_init = [U_init]
    grids = [t_grid, x_grid, 0]
    parameters = [mu]

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
    U_light = np.array(fields_history)[:, 0]
    t_light = time_grid

    pcm = plt.pcolormesh(x_grid, t_light, U_light, cmap='jet', shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$u(x, t)$', rotation=0, size=20, labelpad=-27, y=1.1)
    plt.xlim([x_grid[0], x_grid[-1]])
    plt.xlabel('$x$', size='20')
    plt.ylabel('$t$', size='20')
    #plt.grid(linestyle='--', alpha=0.5)
    plt.savefig('kdV.png', dpi=300)
    plt.close()