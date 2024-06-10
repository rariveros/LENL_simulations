from time_integrators import *

if __name__ == '__main__':

    # Definiendo parámetros
    project_name = '/HO'
    disc = 'C:/'                                        # DISCO DE TRABAJO
    route = 'mnustes_science/simulation_data/FD'        # CARPETA DE TRABAJO
    eq = 'choro_dinger'                                        # ECUACION
    t_rate = 10                                        # CADA CUANTAS ITERACIONES GUARDA
    dt = 0.002
    T = 100
    dx = 0.1 #en milimetros
    alpha = 1  #5.721
    V0 = 1
    phase = np.pi * 0.25

    # Definición de la grilla
    [tmin, tmax, dt] = [0, T, dt]
    [xmin, xmax, dx] = [-10, 10, dx]
    t_grid = np.arange(tmin, tmax + dt, dt)
    x_grid = np.arange(xmin, xmax, dx)
    T = tmax
    Nt = t_grid.shape[0]
    Nx = x_grid.shape[0]

    # Initial Conditions Pattern
    U_1_init = (1 / np.pi) ** 0.25 * np.exp(- x_grid ** 2 / 2) * np.cos(phase)
    U_2_init = (1 / np.pi) ** 0.25 * np.exp(- x_grid ** 2 / 2) * np.sin(phase)

    # Empaquetamiento de parametros, campos y derivadas para integración
    L = xmax - xmin
    D2 = sparse_DD(Nx, dx)
    operators = np.array([D2])
    fields_init = [U_1_init, U_2_init]
    grids = [t_grid, x_grid, 0]
    V = V0 * x_grid ** 2
    parameters = [alpha, V]

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
    arg_light_1 =  (2 * np.pi + arg_light_1) * (arg_light_1 < 0) + arg_light_1 * (arg_light_1 > 0)
    analytical_signal_1 = hilbert(U1_light[-1, :])
    amplitude_envelope_1 = np.abs(analytical_signal_1)

    x = []
    x_sqrt = []
    for i in range(len(t_light)):
        x_i = integrate.simpson(np.conjugate(U_complex[i, :]) * x_grid * U_complex[i, :], x_grid)
        x_i_sqrt = integrate.simpson(np.conjugate(U_complex[i, :]) * x_grid ** 2 * U_complex[i, :], x_grid)
        x.append(x_i)
        x_sqrt.append(x_i_sqrt)
    x = np.array(x)
    x_sqrt = np.array(x_sqrt)
    plt.plot(t_light, x_sqrt - x ** 2)
    plt.show()
    plt.close()

    # Guardando datos
    file = disc + route + project_name
    subfile = "/tset_01"
    parameters_np = np.array([alpha, V0])

    if not os.path.exists(file + subfile):
        os.makedirs(file + subfile)
    np.savetxt(file + subfile + '/field_real.txt', U1_light, delimiter=',')
    np.savetxt(file + subfile + '/field_img.txt', U2_light, delimiter=',')
    np.savetxt(file + subfile + '/parameters.txt', parameters_np, delimiter=',')
    np.savetxt(file + subfile + '/final_envelope.txt', amplitude_envelope_1, delimiter=',')
    np.savetxt(file + subfile + '/X.txt', x_grid, delimiter=',')
    np.savetxt(file + subfile + '/T.txt', t_light, delimiter=',')

    # Guardando gráficos

    pcm = plt.pcolormesh(x_grid, t_light, modulo_light_1, cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$|A|$', rotation=0, size=20, labelpad=-27, y=1.1)
    plt.xlim([x_grid[0], x_grid[-1]])
    plt.xlabel('$x$', size='20')
    plt.ylabel('$t$', size='20')
    plt.grid(linestyle='--', alpha=0.5)
    plt.savefig(file + subfile + '/module_spacetime.png', dpi=300)
    plt.close()

    pcm = plt.pcolormesh(x_grid, t_light, arg_light_1, cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$\\textrm{arg}(A)$', rotation=0, size=20, labelpad=-20, y=1.1)
    plt.xlim([x_grid[0], x_grid[-1]])
    plt.xlabel('$x$', size='20')
    plt.ylabel('$t$', size='20')
    plt.grid(linestyle='--', alpha=0.5)
    plt.savefig(file + subfile + '/arg_spacetime.png', dpi=300)
    plt.close()

    pcm = plt.pcolormesh(x_grid, t_light, U1_light, cmap=parula_map,vmin=-np.amax(U1_light), vmax=np.amax(U1_light), shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$A_R(x, t)$', rotation=0, size=20, labelpad=-27, y=1.1)
    plt.xlim([x_grid[0], x_grid[-1]])
    plt.xlabel('$x$', size='20')
    plt.ylabel('$t$', size='20')
    plt.grid(linestyle='--', alpha=0.5)
    plt.savefig(file + subfile + '/real_spacetime.png', dpi=300)
    plt.close()

    pcm = plt.pcolormesh(x_grid, t_light, U2_light, cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$A_I(x, t)$', rotation=0, size=20, labelpad=-27, y=1.1)
    plt.xlim([x_grid[0], x_grid[-1]])
    plt.xlabel('$x$', size='20')
    plt.ylabel('$t$', size='20')
    plt.grid(linestyle='--', alpha=0.5)
    plt.savefig(file + subfile + '/img_spacetime.png', dpi=300)
    plt.close()