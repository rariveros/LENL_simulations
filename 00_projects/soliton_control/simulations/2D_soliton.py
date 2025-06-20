from back_process import *
from functions import *
from time_integrators import *
from time_integrators import *

if __name__ == '__main__':

    # Definiendo parámetros
    project_name = '/soliton_control/2D'
    disc = 'D:/'
    eq = 'pdnlS_2D'

    alpha = 1 #4 * 6.524
    beta = 1
    gamma_0 = 0.12
    nu = -0.07
    sigma = 15
    mu = 0.1

    gamma_str = str(int(gamma_0 * 1000) * 0.001)
    nu_str = str(int(nu * 1000) * 0.001)
    mu_str = str(int(mu * 1000) * 0.001)
    print('gamma = ' + gamma_str)
    print('nu = ' + nu_str)
    print('mu = ' + mu_str)

    # Definición de la grilla
    Lx = 80
    Ly = 80
    ti = 0
    tf = 1000

    t_rate = 25

    dx = 0.5
    dy = 0.5
    dt = 0.05

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

    c = 0.5
    Z = np.exp((- (X ** 2 + Y ** 2) / (2 * sigma ** 2)) * (1 + 1j * c)) #np.ones((Nx, Ny)) #

    file = disc + 'mnustes_science/simulation_data/FD' + project_name
    subfile = nombre_pndls_gaussian(gamma_0, mu, nu, sigma) #+ '/a=' + A_str
    if not os.path.exists(file + subfile):
        os.makedirs(file + subfile)
    pcm = plt.pcolormesh(x_grid, y_grid, np.abs(Z), cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$\gamma(\\vec{x})$', rotation=0, size=20, labelpad=-27, y=1.15)
    plt.xlabel('$y$', size='20')
    plt.ylabel('$x$', size='20')
    plt.grid(linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(file + subfile + '/forcing.png', dpi=300)
    plt.close()

    # Initial Conditions
    delta = np.sqrt(- nu + np.sqrt(gamma_0 ** 2 - mu ** 2))
    x_0 = 0
    y_0 = 0
    #U_1_init = U_2_init = 0.05 * np.random.rand(Ny, Nx)
    U_1_init = (np.sqrt(2) * delta / np.cosh(delta * (np.sqrt((X - x_0) ** 2 + (Y - y_0) ** 2) / np.sqrt(alpha)))) * np.real(np.exp(-1j * 0.5 * np.arccos(mu / gamma_0)))
    U_2_init = (np.sqrt(2) * delta / np.cosh(delta * (np.sqrt((X - x_0) ** 2 + (Y - y_0) ** 2) / np.sqrt(alpha)))) * np.imag(np.exp(-1j * 0.5 * np.arccos(mu / gamma_0)))

    # Empaquetamiento de parametros, campos y derivadas para integración
    L = xmax - xmin
    D2x = sparse_DD_neumann(Nx, dx)
    D2y = sparse_DD_neumann(Ny, dy)
    operators = np.array([D2x, D2y])
    fields_init = [U_1_init, U_2_init]
    grids = [t_grid, x_grid, y_grid]
    gamma = [np.real(gamma_0 * Z), np.imag(gamma_0 * Z)]

    parameters = [alpha, beta, gamma, mu, nu]

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

    # Definiendo variables finales
    modulo_light_1 = np.absolute(U_complex)
    arg_light_1 = np.angle(U_complex)
    arg_light_1 = (2 * np.pi + arg_light_1) * (arg_light_1 < 0) + arg_light_1 * (arg_light_1 > 0)
    analytical_signal_1 = hilbert(U1_light[-1, :])
    amplitude_envelope_1 = np.abs(analytical_signal_1)

    # Guardando datos
    parameters_np = np.array([alpha, beta, gamma_0, mu, nu])
    U1_light_reshaped = U1_light.reshape(U1_light.shape[0], -1)
    U2_light_reshaped = U2_light.reshape(U2_light.shape[0], -1)

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax))
    ax.grid()
    cax = ax.pcolormesh(x_grid, y_grid, modulo_light_1[0], vmin=0, vmax=0.55,cmap="turbo", shading='auto')
    cbar = fig.colorbar(cax)
    cbar.set_label('$|A(x, y)|$', rotation=0, size=20, labelpad=-27, y=1.13)
    def animate(i):
        cax.set_array(modulo_light_1[i].flatten())
    anim = FuncAnimation(fig, animate, interval=50, frames=len(fields_history) - 1)

    writervideo = animation.FFMpegWriter(fps=60)
    #writervideo = 'imagemagick'
    anim.save(file + subfile + '/soliton_2D.mp4', writer=writervideo, dpi=300)
    plt.close()

