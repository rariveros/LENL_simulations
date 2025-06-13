from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':

    # Definiendo parámetros
    project_name = '/LLG_DMI'
    disc = 'D:/'                                        # DISCO DE TRABAJO
    route = 'mnustes_science/simulation_data/FD'
    eq = 'LLG_DMI'

    t_rate = 2
    A = 4.0
    D = 0.0
    alpha = 0.01

    hx = 0.0
    hy = 0.0
    hz = 1.0

    Kx = 0.0
    Ky = 20.0
    Kz = 0.0

    dist = 50 #30.0
    sigma = 15 #6
    phi = np.pi

    w0 = np.sqrt(hz * (hz + Ky))
    nu = 0.02
    omega = 2 * (w0 + nu) #9.1626

    h = [hx, hy, hz]
    K = [Kx, Ky, Kz]
    L_x = 360
    # Definición de la grilla
    [tmin, tmax, dt] = [0, 2000, 0.02]
    [xmin, xmax, dx] = [- L_x/2, L_x/2, 1.0]
    t_grid = np.arange(tmin, tmax + dt, dt)
    x_grid = np.arange(xmin, xmax, dx)
    T = tmax
    Nt = t_grid.shape[0]
    Nx = x_grid.shape[0]
    dh_0 = 0.37
    dh = dh_0 * np.real((np.exp(- (x_grid - dist / 2) ** 2 / (2 * sigma ** 2)) + np.exp(1j * phi) * np.exp( - (x_grid + dist / 2) ** 2 / (2 * sigma ** 2))))

    mu = alpha * (2 * hz + Ky)
    gamma = dh_0 * (2 * hz + Ky) / omega
    print("nu = " + str(nu))
    print("mu = " + str(mu / 2))
    print("gamma = " + str(gamma / 4))

    # Initial Conditions
    m1_init = m2_init = 0.05 * (np.random.rand(Nx) - 0.5)
    m3_init = np.sqrt(1 - m1_init ** 2 - m2_init ** 2)

    # Empaquetamiento de parametros, campos y derivadas para integración
    L = xmax - xmin
    D1 = sparse_D_neumann(Nx, dx) #sparse_D_neumann(Nx, dx)
    D2 = sparse_DD_neumann(Nx, dx) #sparse_DD_periodic(Nx, dx)
    operators = np.array([D1, D2])
    fields_init = [m1_init, m2_init, m3_init]
    grids = [t_grid, x_grid, 0]
    parameters = [A, D, gamma, alpha, h, dh, omega, K]

    # Midiendo tiempo inicial
    now = datetime.datetime.now()
    print('Hora de Inicio: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
    time_init = time.time()

    final_fields, fields_history, time_grid = RK5_FD(eq, fields_init, parameters, grids, dt, Nt, operators, t_rate)

    now = datetime.datetime.now()
    print('Hora de Término: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
    time_fin = time.time()
    print(str(time_fin - time_init) + ' seg')

    # Reobteniendo campos
    m1_light = np.array(fields_history)[:, 0]
    m2_light = np.array(fields_history)[:, 1]
    m3_light = np.array(fields_history)[:, 2]

    m1_final = final_fields[0]
    m2_final = final_fields[1]
    m3_final = final_fields[2]

    # Guardando datos
    file = disc + route + project_name
    subfile = "/test"
    if not os.path.exists(file + subfile):
        os.makedirs(file + subfile)
    np.savetxt(file + '/m_z.txt', m3_light, delimiter=',')
    np.savetxt(file + '/m_z_final.txt', m3_final, delimiter=',')
    np.savetxt(file + '/x_grid.txt', x_grid, delimiter=',')
    np.savetxt(file + '/t_grid.txt', time_grid, delimiter=',')

    # Graficando
    pcm = plt.pcolormesh(x_grid, time_grid, m1_light, cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$m_x(x,t)$', rotation=0, size=25, labelpad=-27, y=1.12)
    plt.xlim([x_grid[0], x_grid[-1]])
    plt.ylim([1900, t_grid[-1]])
    plt.xlabel('$x$', size='25')
    plt.ylabel('$t$', size='25')
    cbar.ax.tick_params(labelsize=15)
    plt.grid(linestyle='--', alpha=0.25)
    plt.savefig(file + '/mx_final.png', dpi=300)
    plt.tight_layout()
    plt.close()

    pcm = plt.pcolormesh(x_grid, time_grid, m2_light, cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$m_y(x,t)$', rotation=0, size=25, labelpad=-27, y=1.12)
    plt.xlim([x_grid[0], x_grid[-1]])
    plt.xlabel('$x$', size='25')
    plt.ylabel('$t$', size='25')
    cbar.ax.tick_params(labelsize=15)
    plt.grid(linestyle='--', alpha=0.25)
    plt.savefig(file + '/my_final.png', dpi=300)
    plt.tight_layout()
    plt.close()

    pcm = plt.pcolormesh(x_grid, time_grid, m3_light, cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$m_z(x,t)$', rotation=0, size=25, labelpad=-27, y=1.12)
    plt.xlim([x_grid[0], x_grid[-1]])
    plt.xlabel('$x$', size='25')
    plt.ylabel('$t$', size='25')
    cbar.ax.tick_params(labelsize=15)
    plt.grid(linestyle='--', alpha=0.25)
    plt.savefig(file + '/mz_final.png', dpi=300)
    plt.tight_layout()
    plt.close()

    pcm = plt.pcolormesh(x_grid, time_grid, m1_light ** 2 + m2_light ** 2 + m3_light ** 2, cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$|m(x,t)|^{2}$', rotation=0, size=25, labelpad=-27, y=1.12)
    plt.xlim([x_grid[0], x_grid[-1]])
    plt.xlabel('$x$', size='25')
    plt.ylabel('$t$', size='25')
    cbar.ax.tick_params(labelsize=15)
    plt.grid(linestyle='--', alpha=0.25)
    plt.savefig(file + '/m_mod.png', dpi=300)
    plt.tight_layout()
    plt.close()

    fig, ax01 = plt.subplots(1, 1, figsize=(8, 2))
    ax01.plot(x_grid, m1_light[-1, :], color="r", label="$m_x$", lw=3)
    ax01.plot(x_grid, m2_light[-1, :], color="g", label="$m_y$", lw=3)
    ax01.plot(x_grid, m3_light[-1, :], color="b", label="$m_z$", lw=3)
    ax01.set_xlim([x_grid[0], x_grid[-1]])
    ax01.tick_params(labelsize=20, labelleft=True, labelright=False, left=True, right=False)
    ax01.legend(loc="upper right", fontsize=20)
    plt.grid(linestyle='--', alpha=0.5)
    plt.savefig(file + '/pretty_init_profiles.png', dpi=200, bbox_inches='tight')
    plt.close()