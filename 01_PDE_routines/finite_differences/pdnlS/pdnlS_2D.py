from back_process import *
from functions import *
from time_integrators import *


from functions import *
from back_process import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':

    # Definiendo parámetros
    project_name = '/pdnlS_2D'
    disc = 'C:/'
    eq = 'pdnlS_2D'

    alpha = 1
    beta = 1
    gamma_0 = 0.15
    nu = 0.4
    sigma = 10
    mu = 0.1

    gamma_str = str(int(gamma_0 * 1000) * 0.001)
    nu_str = str(int(nu * 1000) * 0.001)
    mu_str = str(int(mu * 1000) * 0.001)
    print('gamma = ' + gamma_str)
    print('nu = ' + nu_str)
    print('mu = ' + mu_str)

    # Definición de la grilla
    Lx = 100
    Ly = 100
    ti = 0
    tf = 250

    dx = 0.25
    dy = 0.25
    dt = 0.02

    [tmin, tmax, dt] = [0, tf, dt]
    [xmin, xmax, dx] = [- Lx / 2,  Lx / 2, dx]
    [ymin, ymax, dy] = [- Ly / 2, Ly / 2, dy]
    t_grid = np.arange(tmin, tmax + dt, dt)
    x_grid = np.arange(xmin, xmax, dx)
    y_grid = np.arange(ymin, ymax, dy)
    T = tmax
    Nt = t_grid.shape[0]
    Nx = x_grid.shape[0]
    Ny = y_grid.shape[0]
    X, Y = np.meshgrid(x_grid, y_grid)


    # Initial Conditions
    U_1_init = U_2_init = 0.01 * np.random.rand(Ny, Nx)

    # Empaquetamiento de parametros, campos y derivadas para integración
    L = xmax - xmin
    D2x = sparse_DD_neumann(Nx, dx)
    D2y = sparse_DD_neumann(Ny, dy)
    operators = np.array([D2x, D2y])
    fields_init = [U_1_init, U_2_init]
    grids = [t_grid, x_grid, y_grid]
    gamma = [gamma_0, 0]

    parameters = [alpha, beta, gamma, mu, nu]

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
    subfile = nombre_pndls_gaussian(gamma_0, mu, nu, sigma) #+ '/a=' + A_str
    parameters_np = np.array([alpha, beta, gamma_0, mu, nu])
    if not os.path.exists(file + subfile):
        os.makedirs(file + subfile)
    U1_light_reshaped = U1_light.reshape(U1_light.shape[0], -1)
    U2_light_reshaped = U2_light.reshape(U2_light.shape[0], -1)
    np.savetxt(file + subfile + '/field_real.txt', U1_light_reshaped, delimiter=',')
    np.savetxt(file + subfile + '/field_img.txt', U2_light_reshaped, delimiter=',')
    np.savetxt(file + subfile + '/parameters.txt', parameters_np, delimiter=',')
    np.savetxt(file + subfile + '/final_envelope.txt', amplitude_envelope_1, delimiter=',')
    np.savetxt(file + subfile + '/X.txt', x_grid, delimiter=',')
    np.savetxt(file + subfile + '/T.txt', t_light, delimiter=',')

    nu_positive_grid = np.arange(0, 2, 0.01)
    nu_negative_grid = - np.flip(nu_positive_grid)
    nu_grid = np.append(nu_negative_grid, nu_positive_grid)
    plt.plot(nu_positive_grid, np.sqrt(nu_positive_grid ** 2 + mu ** 2), c='k', linestyle='--')
    plt.fill_between(nu_positive_grid, np.ones(len(nu_positive_grid)) * mu,
                     np.sqrt(nu_positive_grid ** 2 + mu ** 2),
                     facecolor=(92 / 255, 43 / 255, 228 / 255, 0.4))
    plt.plot(nu_negative_grid, np.sqrt(nu_negative_grid ** 2 + mu ** 2), c='k', linestyle='--')
    plt.fill_between(nu_negative_grid, np.ones(len(nu_negative_grid)) * mu,
                     np.sqrt(nu_negative_grid ** 2 + mu ** 2),
                     facecolor=(0, 1, 0, 0.4))
    plt.plot(nu_grid, np.ones(len(nu_grid)) * mu, c='k', linestyle='--')
    plt.fill_between(nu_grid, 2, np.sqrt(nu_grid ** 2 + mu ** 2),
                     facecolor=(1, 0, 0, 0.4))
    plt.fill_between(nu_grid, np.ones(len(nu_grid)) * mu, 0,
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

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.set(xlim=(x_grid[0], x_grid[-1]), ylim=(y_grid[0], y_grid[-1]))
    cax = ax.pcolormesh(x_grid, y_grid, fields_history[0][0], vmin=np.amin(final_fields[0]), vmax=np.amax(final_fields[0]), cmap=parula_map, shading='auto')
    cbar = fig.colorbar(cax)
    cbar.set_label('$A_R(x, y)$', rotation=0, size=20, labelpad=-27, y=1.13)
    def animate(i):
        cax.set_array(fields_history[i][0].flatten())
    anim = FuncAnimation(fig, animate, interval=50, frames=len(fields_history) - 1)

    writervideo = animation.FFMpegWriter(fps=60)
    #writervideo = 'imagemagick'
    anim.save(file + subfile + '/pdnlS_2D.mp4', writer=writervideo, dpi=300)
    plt.close()

    pcm = plt.pcolormesh(x_grid, y_grid, final_fields[0], vmin=np.amin(final_fields[0]), vmax=np.amax(final_fields[0]), cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$A_R(x, y)$', rotation=0, size=20, labelpad=-27, y=1.15)
    plt.xlabel('$y$', size='20')
    plt.ylabel('$x$', size='20')
    plt.grid(linestyle='--', alpha=0.5)
    plt.savefig(file + subfile + '/final_profile.png', dpi=300)
    plt.close()
