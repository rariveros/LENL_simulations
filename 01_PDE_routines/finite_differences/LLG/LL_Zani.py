from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':

    # Definiendo parámetros
    project_name = '/LL_Zani_homo'
    disc = 'C:/'
    eq = 'LL_Zani'

    alpha = 0.03
    C_ani = 0.2
    sigma = 70

    L_y = L_x = 500

    # Definición de la grilla
    [tmin, tmax, dt] = [0, 600, 0.075]
    [xmin, xmax, dx] = [- L_x/2, L_x/2, 1]
    [ymin, ymax, dy] = [- L_y/2, L_y/2, 1]
    t_grid = np.arange(tmin, tmax + dt, dt)
    x_grid = np.arange(xmin, xmax, dx)
    y_grid = np.arange(ymin, ymax, dy)
    X, Y = np.meshgrid(x_grid, y_grid)
    T = tmax
    Nt = t_grid.shape[0]
    Nx = x_grid.shape[0]
    Ny = y_grid.shape[0]

    # Initial Conditions
    m2_init = 0.01 * np.random.rand(Nx, Ny)
    m3_init = 0.01 * np.random.rand(Nx, Ny)
    m1_init = 1 - (m2_init) ** 2 - (m3_init) ** 2

    # Empaquetamiento de parametros, campos y derivadas para integración
    L = xmax - xmin
    D2 = sparse_DD(Nx, dx)
    operators = np.array([D2])
    fields_init = [m1_init, m2_init, m3_init]
    grids = [t_grid, x_grid, y_grid]
    parameters = [alpha, C_ani] #* np.exp(- (X ** 2 + Y ** 2) / (2 * sigma ** 2))

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
    m1_light = np.array(fields_history)[:, 0]
    m2_light = np.array(fields_history)[:, 1]
    m3_light = np.array(fields_history)[:, 2]


    m1_final = final_fields[0]
    m2_final = final_fields[1]
    m3_final = final_fields[2]

    # Guardando datos
    file = disc + 'Users/mnustes_science/PT_fluids/mnustes_science/simulation_data/FD' + project_name
    if not os.path.exists(file):
        os.makedirs(file)
    m3_light_reshaped = m3_light.reshape(m3_light.shape[0], -1)
    #loaded_arr = np.loadtxt("geekfile.txt")
    #load_original_arr = loaded_arr.reshape(loaded_arr.shape[0], loaded_arr.shape[1] // arr.shape[2], arr.shape[2])
    np.savetxt(file + '/m_z.txt', m3_light_reshaped, delimiter=',')
    np.savetxt(file + '/m_z_final.txt', m3_final, delimiter=',')
    np.savetxt(file + '/x_grid.txt', x_grid, delimiter=',')
    np.savetxt(file + '/y_grid.txt', y_grid, delimiter=',')
    np.savetxt(file + '/t_grid.txt', time_grid, delimiter=',')

    # Graficando
    pcm = plt.pcolormesh(x_grid, y_grid, m3_final, cmap='plasma', shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$m_z(x,y)$', rotation=0, size=25, labelpad=-27, y=1.1)
    plt.xlim([x_grid[0], x_grid[-1]])
    #plt.xticks([-200, -100, 0, 100, 100], fontsize=20)
    #plt.yticks([-200, -100, 0, 100, 100], fontsize=20)

    plt.xlabel('$x$', size='25')
    plt.ylabel('$y$', size='25')

    #plt.ylim([-200, 200])
    cbar.ax.tick_params(labelsize=15)
    plt.grid(linestyle='--', alpha=0.25)
    plt.savefig(file + '/mz_final.png', dpi=300)
    plt.tight_layout()
    plt.close()

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.set(xlim=(-L_x / 2, L_x / 2), ylim=(-L_y / 2, L_y / 2))
    cax = ax.pcolormesh(x_grid, y_grid, fields_history[0][2], vmin=-1, vmax=1, cmap='plasma', shading='auto')
    cbar = fig.colorbar(cax)
    cbar.set_label('$m_z(x, y)$', rotation=0, size=20, labelpad=-27, y=1.13)
    def animate(i):
        cax.set_array(fields_history[i][2].flatten())
    anim = FuncAnimation(
        fig, animate, interval=10, frames=len(fields_history) - 1)

    FFwriter = animation.FFMpegWriter()
    anim.save(file + '/mzdynamics.gif', writer='imagemagick', fps=60, dpi=300)
    plt.close()
