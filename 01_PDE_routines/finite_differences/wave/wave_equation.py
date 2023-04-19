from back_process import *
from back_process import *
from functions import *
from time_integrators import *

if __name__ == '__main__':
    eq = 'wave'
    L_x = 20
    L_y = L_x
    dx = 0.02
    dy = dx
    x_grid = np.arange(-L_x/2, L_x/2, dx)
    y_grid = np.arange(-L_y / 2, L_y / 2, dy)
    Nx = len(x_grid)
    Ny = len(y_grid)
    X, Y = np.meshgrid(x_grid, y_grid)
    n = 5
    m = 3
    c = 1
    sigma = 0.5
    U_1_init = np.exp(-(X ** 2 + Y ** 2) / (2 * sigma ** 2))#np.cos(2 * n * np.pi * X / L_x) * np.cos(2 * m * np.pi * Y / L_y)
    U_2_init = np.full_like(U_1_init, 0)

    dt = 0.01
    T = 10
    t_grid = np.arange(0, T, dt)
    Nt = len(t_grid)
    D2 = sparse_DD_neumann(Nx, dx)
    operators = np.array([D2])
    parameters = [c]
    fields_init = [U_1_init, U_2_init]
    grids = [t_grid, x_grid, y_grid]

    # Midiendo tiempo inicial
    now = datetime.datetime.now()
    print('Hora de Inicio: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
    time_init = time.time()

    final_fields, fields_history, time_grid = RK4_FD(eq, fields_init, parameters, grids, dt, Nt, operators)

    now = datetime.datetime.now()
    print('Hora de Término: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
    time_fin = time.time()
    print(str(time_fin - time_init) + ' seg')

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.set(xlim=(-L_x / 2, L_x / 2), ylim=(-L_y / 2, L_y / 2))
    cax = ax.pcolormesh(x_grid, y_grid, fields_history[0][0], vmin=-1, vmax=1, cmap=parula_map, shading='auto')
    cbar = fig.colorbar(cax)
    cbar.set_label('$u(x, y)$', rotation=0, size=20, labelpad=-27, y=1.13)
    def animate(i):
        cax.set_array(fields_history[i][0].flatten())
    anim = FuncAnimation(
        fig, animate, interval=10, frames=len(fields_history) - 1)

    FFwriter = animation.FFMpegWriter()
    anim.save('wave.gif', writer='imagemagick', fps=60, dpi=300)
    plt.close()

    pcm = plt.pcolormesh(x_grid, y_grid, final_fields[0], vmin=-1, vmax=1, cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$u(x, y)$', rotation=0, size=20, labelpad=-27, y=1.15)
    plt.xlabel('$y$', size='20')
    plt.ylabel('$x$', size='20')
    plt.grid(linestyle='--', alpha=0.5)
    plt.savefig('final_profile.png', dpi=300)
    plt.close()
