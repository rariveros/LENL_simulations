from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':

    # Definiendo parámetros
    project_name = '/LLGZ_01'
    disc = 'C:/'
    eq = 'LLGZ_01'

    alpha = 0.1
    h_1 = - 0.3
    sigma = 3
    A = 1
    g = -0.05
    L_x = 60
    h_d_0 = 0.25
    # Definición de la grilla
    [tmin, tmax, dt] = [0, 250, 0.0005]
    [xmin, xmax, dx] = [- L_x/2, L_x/2, 0.4]
    t_grid = np.arange(tmin, tmax + dt, dt)
    x_grid = np.arange(xmin, xmax, dx)
    T = tmax
    Nt = t_grid.shape[0]
    Nx = x_grid.shape[0]

    h_d =h_d_0 #* (np.exp(-(x_grid + 10) ** 2 / (2 * sigma ** 2)) - np.exp(-(x_grid - 10) ** 2 / (2 * sigma ** 2)))

    print("mu = " + str(-g))
    print("nu = " + str(-h_1 - h_d_0 / 2))
    print("gamma = " + str(h_d / 2))
    # Initial Conditions
    m1_init = m2_init = 0.01 * np.random.rand(Nx)
    m3_init = np.sqrt(1 - m1_init ** 2 - m2_init ** 2)

    #plt.plot(x_grid, m1_init, label="$m_x$")
    #plt.plot(x_grid, m2_init, label="$m_y$")
    #plt.plot(x_grid, m3_init, label="$m_z$")
    #plt.legend()
    #plt.show()

    # Empaquetamiento de parametros, campos y derivadas para integración
    L = xmax - xmin
    D2 = sparse_DD(Nx, dx)
    operators = np.array([D2])
    fields_init = [m1_init, m2_init, m3_init]
    grids = [t_grid, x_grid, 0]
    parameters = [alpha, g, h_1, h_d, A]

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
    np.savetxt(file + '/m_z.txt', m3_light, delimiter=',')
    np.savetxt(file + '/m_z_final.txt', m3_final, delimiter=',')
    np.savetxt(file + '/x_grid.txt', x_grid, delimiter=',')
    np.savetxt(file + '/t_grid.txt', time_grid, delimiter=',')

    # Graficando
    pcm = plt.pcolormesh(x_grid, time_grid, m1_light, cmap='plasma', shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$m_x(x,t)$', rotation=0, size=25, labelpad=-27, y=1.12)
    plt.xlim([x_grid[0], x_grid[-1]])
    plt.xlabel('$x$', size='25')
    plt.ylabel('$t$', size='25')
    cbar.ax.tick_params(labelsize=15)
    plt.grid(linestyle='--', alpha=0.25)
    plt.savefig(file + '/mx_final.png', dpi=300)
    plt.tight_layout()
    plt.close()

    pcm = plt.pcolormesh(x_grid, time_grid, m2_light, cmap='plasma', shading='auto')
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

    pcm = plt.pcolormesh(x_grid, time_grid, m3_light, cmap='plasma', shading='auto')
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

    pcm = plt.pcolormesh(x_grid, time_grid, m1_light ** 2 + m2_light ** 2 + m3_light ** 2, cmap='jet', shading='auto')
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
    #fig, ax = plt.subplots(figsize=(5, 4))
    #ax.set(xlim=(-L_x / 2, L_x / 2), ylim=(-L_y / 2, L_y / 2))
    #cax = ax.pcolormesh(x_grid, y_grid, fields_history[0][2], vmin=-1, vmax=1, cmap='plasma', shading='auto')
    #cbar = fig.colorbar(cax)
    #cbar.set_label('$m_z(x, y)$', rotation=0, size=20, labelpad=-27, y=1.13)
    #def animate(i):
    #    cax.set_array(fields_history[i][2].flatten())
    #anim = FuncAnimation(
    #    fig, animate, interval=10, frames=len(fields_history) - 1)
    #
    #FFwriter = animation.FFMpegWriter()
    #anim.save(file + '/mzdynamics.gif', writer='imagemagick', fps=60, dpi=300)
    #plt.close()
