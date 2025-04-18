import matplotlib.pyplot as plt
from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':

    # Definiendo parámetros
    project_name = '/pdnlS_dimer_tranformed'
    disc = 'C:/'
    route = 'mnustes_science/simulation_data/FD'
    eq = 'pdnlS_dimer_tranformed'
    t_rate = 1

    nu = 0.0
    mu = 0.1
    Gammas = [0.13] #np.arange(0.09, 0.2, 0.02)
    K = [0.1]#np.arange(0.00, 0.51, 0.02)
    g = 1.0

    # Definición de la grilla
    [tmin, tmax, dt] = [0, 1000, 0.2]
    [xmin, xmax, dx] = [0, 1, 1]
    t_grid = np.arange(tmin, tmax + dt, dt)
    x_grid = np.arange(xmin, xmax, dx)
    T = tmax
    Nt = t_grid.shape[0]
    Nx = x_grid.shape[0]
    operators = [0]
    ti = 0.
    U_init = 0.01 * (np.random.rand(Nx) + 1j * np.random.rand(Nx))
    V_init = 0.0 * (np.random.rand(Nx) + 1j * np.random.rand(Nx))
    for gamma in Gammas:
        for k in K:
            Delta_str = f"{nu:.{4}f}"
            gamma_str = f"{mu:.{4}f}"
            Omega_str = f"{gamma:.{4}f}"
            k_str = f"{k:.{4}f}"

            print("####### " + k_str + " #######")

            # Empaquetamiento de parametros, campos y derivadas para integración
            fields_init = [U_init, V_init]
            grids = [t_grid, x_grid, 0]

            parameters_np = np.array([nu, mu, gamma, k, g])

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
            subfile = "/nu=" + Delta_str + "/mu=" + gamma_str + "/k=" + k_str + "/gamma=" + Omega_str
            if not os.path.exists(file + subfile):
                os.makedirs(file + subfile)

            np.savetxt(file + subfile + '/U1.txt', U, delimiter=',')
            np.savetxt(file + subfile + '/V1.txt', V, delimiter=',')
            np.savetxt(file + subfile + '/parameters.txt', parameters_np, delimiter=',')
            np.savetxt(file + subfile + '/T.txt', t_grid, delimiter=',')

            lightness = 1
            U_light = U[0::lightness]
            V_light = V[0::lightness]
            t_light = time_grid[0::lightness]

            fig, (ax1, ax2) = plt.subplots(2)
            ax1.plot(t_light, np.real(U_light), c="b")
            ax1.plot(t_light, np.imag(U_light), c="r")
            plt.xlabel('$x$', size='20')
            plt.ylabel('$t$', size='20')
            ax1.set_xlim([ti * t_light[-1], t_light[-1]])
            ax1.grid(linestyle='--', alpha=0.5)

            ax2.plot(t_light, np.real(V_light), c="b")
            ax2.plot(t_light, np.imag(V_light), c="r")
            plt.xlabel('$x$', size='20')
            plt.ylabel('$t$', size='20')
            ax2.set_xlim([ti * t_light[-1], t_light[-1]])
            ax2.grid(linestyle='--', alpha=0.5)
            plt.savefig(file + subfile + "/timeseries.png", dpi=300)
            plt.close()

            U_init = U_light[-1] #+ 0.01 * (np.random.rand(Nx) + 1j * np.random.rand(Nx))
            V_init = V_light[-1] #+ 0.01 * (np.random.rand(Nx) + 1j * np.random.rand(Nx))