import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':

    # Definiendo parámetros
    project_name = '/mathieu/single'
    disc = 'C:/'
    route = 'mnustes_science/simulation_data/FD'
    eq = 'mathieu_single'
    t_rate = 100

    w_0 = 1
    w_i = 0.25
    gamma_i = 0.05
    delta = 0

    ti = 0.

    sigma = 2
    # Definición de la grilla
    [tmin, tmax, dt] = [0, 1000, 0.005]
    [xmin, xmax, dx] = [0, 1, 1]
    t_grid = np.arange(tmin, tmax + dt, dt)
    x_grid = np.arange(xmin, xmax, dx)
    T = tmax
    Nt = t_grid.shape[0]
    Nx = x_grid.shape[0]


    # Initial Conditions Pattern
    U_init = 0.0001 * np.random.rand(Nx) #np.exp(-x_grid ** 2 / (2 * sigma ** 2) ** 2)
    V_init = 0 * np.random.rand(Nx)

    #DD = sparse_DD_neumann(Nx, dx)
    #operators = np.array([DD])
    operators = [0]

    #interaction = np.array([[0, 1], [1, 0]])
    #operators =[interaction]

    # Empaquetamiento de parametros, campos y derivadas para integración
    fields_init = [U_init, V_init]
    grids = [t_grid, x_grid, 0]

    parameters_np = np.array([w_0, w_i, gamma_i, delta])

    # Midiendo tiempo inicial
    now = datetime.datetime.now()
    print('Hora de Inicio: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
    time_init = time.time()

    final_fields, fields_history, time_grid = RK4_FD(eq, fields_init, parameters_np, grids, dt, Nt, operators, t_rate)

    now = datetime.datetime.now()
    print('Hora de Término: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
    time_fin = time.time()
    print(str(time_fin - time_init) + ' seg')

    # Reobteniendo campos
    U = np.array(fields_history)[:, 0]
    V = np.array(fields_history)[:, 1]

    # Guardando datos
    file = disc + route + project_name
    subfile = "/test"
    if not os.path.exists(file + subfile):
        os.makedirs(file + subfile)

    np.savetxt(file + subfile + '/U.txt', U, delimiter=',')
    np.savetxt(file + subfile + '/V.txt', V, delimiter=',')
    np.savetxt(file + subfile + '/parameters.txt', parameters_np, delimiter=',')
    np.savetxt(file + subfile + '/T.txt', t_grid, delimiter=',')

    lightness = 1
    U_light = U[0::lightness]
    V_light = V[0::lightness]
    t_light = time_grid[0::lightness]
    forcing = (gamma_i / 2) * (np.cos(w_i * np.array(t_light)) + np.cos((1 + delta) * w_i * np.array(t_light)))

    fig, (ax1, ax2) = plt.subplots(2)
    ax1.plot(t_light, U_light, c="k")
    plt.xlabel('$x$', size='20')
    plt.ylabel('$t$', size='20')
    ax1.set_xlim([ti * t_light[-1], t_light[-1]])
    ax1.grid(linestyle='--', alpha=0.5)

    ax2.plot(t_light, forcing, c="k")
    plt.xlabel('$x$', size='20')
    plt.ylabel('$t$', size='20')
    ax2.set_xlim([ti * t_light[-1], t_light[-1]])
    ax2.grid(linestyle='--', alpha=0.5)
    plt.show()
    plt.close()





