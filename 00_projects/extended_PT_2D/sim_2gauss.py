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
    gamma_0 = 0.28
    nu = 0.32
    sigma = 4
    mu = 0.1

    gamma_str = str(int(gamma_0 * 1000) * 0.001)
    nu_str = str(int(nu * 1000) * 0.001)
    mu_str = str(int(mu * 1000) * 0.001)
    print('gamma = ' + gamma_str)
    print('nu = ' + nu_str)
    print('mu = ' + mu_str)

    # Definición de la grilla
    Lx = 120
    Ly = 120
    ti = 0
    tf = 2500

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
    R = 15
    N = 6
    Zs = []
    for i in range(N):
        x_i = R * np.cos(2 * np.pi * i / N)
        y_i = R * np.sin(2 * np.pi * i / N)
        Z_i = ((-1) ** i) * np.exp(-((X - x_i) ** 2 + (Y - y_i) ** 2) / (2 * sigma ** 2))
        Zs.append(Z_i)
    Zs = np.array(Zs)
    Z = Zs[0]
    for i in range(1, len(Zs)):
        Z = Z + Zs[i]

    file = disc + 'mnustes_science/simulation_data/FD' + project_name
    subfile = nombre_pndls_gaussian(gamma_0, mu, nu, sigma) #+ '/a=' + A_str
    if not os.path.exists(file + subfile):
        os.makedirs(file + subfile)
    pcm = plt.pcolormesh(x_grid, y_grid, Z, cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$\gamma(\\vec{x})$', rotation=0, size=20, labelpad=-27, y=1.15)
    plt.xlabel('$y$', size='20')
    plt.ylabel('$x$', size='20')
    plt.grid(linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(file + subfile + '/forcing.png', dpi=300)
    plt.close()

    # Initial Conditions
    U_1_init = U_2_init = 0.01 * np.random.rand(Ny, Nx)

    # Empaquetamiento de parametros, campos y derivadas para integración
    L = xmax - xmin
    D2x = sparse_DD_neumann(Nx, dx)
    D2y = sparse_DD_neumann(Ny, dy)
    operators = np.array([D2x, D2y])
    fields_init = [U_1_init, U_2_init]
    grids = [t_grid, x_grid, y_grid]
    gamma = [gamma_0 * Z, 0]

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
    ax.set(xlim=(-40, 40), ylim=(-40, 40))
    cax = ax.pcolormesh(x_grid, y_grid, modulo_light_1[0], vmin=0, vmax=0.55,cmap="turbo", shading='auto')
    cbar = fig.colorbar(cax)
    cbar.set_label('$|A(x, y)|$', rotation=0, size=20, labelpad=-27, y=1.13)
    def animate(i):
        cax.set_array(modulo_light_1[i].flatten())
    anim = FuncAnimation(fig, animate, interval=50, frames=len(fields_history) - 1)

    writervideo = animation.FFMpegWriter(fps=60)
    #writervideo = 'imagemagick'
    anim.save(file + subfile + '/pdnlS_2D.mp4', writer=writervideo, dpi=300)
    plt.close()

