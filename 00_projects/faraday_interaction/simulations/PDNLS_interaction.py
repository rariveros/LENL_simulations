from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':

    # Definiendo parámetros
    project_name = '/PDNLS_interaction'
    disc = 'C:/'
    eq = 'PDNLS_interaction'
    t_rate = 100
    d_nu = 0.05
    d_gammas = 0.05


    alpha_1 = 1
    beta_1 = 1
    gamma1_0 = 0.22
    mu_1 = 0.1
    nu_1 = 0.25
    sigma_1 = 8

    alpha_2 = 1
    beta_2 = 1
    gamma2_0 = 0.22
    mu_2 = 0.1
    nu_2 = nu_1 / 5
    sigma_2 = 8

    gamma_str = str(int(gamma1_0 * 1000) * 0.001)
    nu_str = str(int(nu_1 * 1000) * 0.001)
    mu_str = str(int(mu_1 * 1000) * 0.001)
    print('gamma = ' + gamma_str)
    print('nu = ' + nu_str)
    print('mu = ' + mu_str)

    # Definición de la grilla
    [tmin, tmax, dt] = [0, 1500, 0.02]
    [xmin, xmax, dx] = [-50, 50, 0.2]
    t_grid = np.arange(tmin, tmax + dt, dt)
    x_grid = np.arange(xmin, xmax, dx)
    T = tmax
    Nt = t_grid.shape[0]
    Nx = x_grid.shape[0]

    # Initial Conditions Pattern
    U_1_init = 0.01 * np.random.rand(Nx)
    U_2_init = 0.01 * np.random.rand(Nx)
    V_1_init = 0.01 * np.random.rand(Nx)
    V_2_init = 0.01 * np.random.rand(Nx)

    # Empaquetamiento de parametros, campos y derivadas para integración
    L = xmax - xmin
    D2 = sparse_DD_neumann(Nx, dx)
    operators = np.array([D2])
    fields_init = [U_1_init, U_2_init, V_1_init, V_2_init]
    grids = [t_grid, x_grid, 0]
    gamma1_real = gamma1_0 * np.exp(- x_grid ** 2 / (2 * sigma_1 ** 2))
    gamma1_img = 0
    gamma1 = [gamma1_real, gamma1_img]

    gamma2_real = gamma2_0 * np.exp(- x_grid ** 2 / (2 * sigma_2 ** 2))
    gamma2_img = 0
    gamma2 = [gamma2_real, gamma2_img]

    parameters = [alpha_1, beta_1, gamma1, mu_1, nu_1, alpha_2, beta_2, gamma2, mu_2, nu_2]

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
    U2_light = np.array(fields_history)[:, 1]
    U_complex = U1_light + 1j * U2_light
    t_light = time_grid

    V1_light = np.array(fields_history)[:, 2]
    V2_light = np.array(fields_history)[:, 3]
    V_complex = V1_light + 1j * V2_light

    pcm = plt.pcolormesh(x_grid, t_light, U1_light, cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$A_R(x, t)$', rotation=0, size=20, labelpad=-27, y=1.1)
    plt.xlim([x_grid[0], x_grid[-1]])
    plt.xlabel('$x$', size='20')
    plt.ylabel('$t$', size='20')
    plt.grid(linestyle='--', alpha=0.5)
    plt.savefig('real_spacetime_01.png', dpi=300)
    plt.close()

    pcm = plt.pcolormesh(x_grid, t_light, V1_light, cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$B_R(x, t)$', rotation=0, size=20, labelpad=-27, y=1.1)
    plt.xlim([x_grid[0], x_grid[-1]])
    plt.xlabel('$x$', size='20')
    plt.ylabel('$t$', size='20')
    plt.grid(linestyle='--', alpha=0.5)
    plt.savefig('real_spacetime_02.png', dpi=300)
    plt.close()

    pcm = plt.pcolormesh(x_grid, t_light, U1_light + V1_light, cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$A_R + B_R$', rotation=0, size=20, labelpad=-27, y=1.1)
    plt.xlim([x_grid[0], x_grid[-1]])
    plt.xlabel('$x$', size='20')
    plt.ylabel('$t$', size='20')
    plt.grid(linestyle='--', alpha=0.5)
    plt.savefig('real_spacetime_joint.png', dpi=300)
    plt.close()