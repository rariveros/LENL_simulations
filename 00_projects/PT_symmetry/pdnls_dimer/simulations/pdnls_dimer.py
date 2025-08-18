import matplotlib.pyplot as plt
from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':

    # Definiendo parámetros
    project_name = '/pdnlS_dimer_test'
    disc = 'D:/'
    route = 'mnustes_science/simulation_data/FD'
    eq = 'pdnlS_dimer'
    t_rate = 1

    nu = 0.1
    mu = 0.1
    Gammas = np.arange(0.09, 0.2, 0.02)
    K = [0.1]#np.arange(0.00, 0.51, 0.02)
    g = 1.0

    # Definición de la grilla
    [tmin, tmax, dt] = [0, 1000, 0.1]
    [xmin, xmax, dx] = [0, 1, 1]
    t_grid = np.arange(tmin, tmax + dt, dt)
    x_grid = np.arange(xmin, xmax, dx)
    T = tmax
    Nt = t_grid.shape[0]
    Nx = x_grid.shape[0]
    operators = [0]
    ti = 0.
    U1_init = 0.01 * (np.random.rand(Nx) + 1j * np.random.rand(Nx))
    U2_init = np.conjugate(U1_init)
    V1_init = 0.01 * (np.random.rand(Nx) + 1j * np.random.rand(Nx))
    V2_init = np.conjugate(V1_init)
    for gamma in Gammas:
        for k in K:
            nu_str = f"{nu:.{4}f}"
            mu_str = f"{mu:.{4}f}"
            gamma_str = f"{gamma:.{4}f}"
            k_str = f"{k:.{4}f}"

            print("####### " + k_str + " #######")

            # Empaquetamiento de parametros, campos y derivadas para integración
            fields_init = [U1_init, U2_init, V1_init, V2_init]
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
            U1 = np.array(fields_history)[:, 0]
            U2 = np.array(fields_history)[:, 1]
            V1 = np.array(fields_history)[:, 2]
            V2 = np.array(fields_history)[:, 3]

            # Guardando datos
            file = disc + route + project_name
            subfile = "/nu=" + nu_str + "/mu=" + mu_str + "/gamma=" + gamma_str + "/k=" + k_str
            if not os.path.exists(file + subfile):
                os.makedirs(file + subfile)

            np.savetxt(file + subfile + '/U1.txt', U1, delimiter=',')
            np.savetxt(file + subfile + '/U2.txt', U1, delimiter=',')
            np.savetxt(file + subfile + '/V1.txt', V1, delimiter=',')
            np.savetxt(file + subfile + '/V2.txt', V2, delimiter=',')
            np.savetxt(file + subfile + '/parameters.txt', parameters_np, delimiter=',')
            np.savetxt(file + subfile + '/T.txt', t_grid, delimiter=',')

            lightness = 1
            U1_light = U1[0::lightness]
            U2_light = U2[0::lightness]
            V1_light = V1[0::lightness]
            V2_light = V2[0::lightness]
            t_light = time_grid[0::lightness]

            fig, (ax1, ax2) = plt.subplots(2)
            ax1.plot(t_light, np.real(U1_light), c="b")
            ax1.plot(t_light, np.imag(U1_light), c="r")
            plt.xlabel('$x$', size='20')
            plt.ylabel('$t$', size='20')
            ax1.set_xlim([ti * t_light[-1], t_light[-1]])
            ax1.grid(linestyle='--', alpha=0.5)

            ax2.plot(t_light, np.real(V1_light), c="b")
            ax2.plot(t_light, -np.imag(V1_light), c="r")
            plt.xlabel('$x$', size='20')
            plt.ylabel('$t$', size='20')
            ax2.set_xlim([ti * t_light[-1], t_light[-1]])
            ax2.grid(linestyle='--', alpha=0.5)
            plt.savefig(file + subfile + "/timeseries.png", dpi=300)
            plt.close()

            U1_init = U1_light[-1] + 0.01 * (np.random.rand(Nx) + 1j * np.random.rand(Nx))
            U2_init = np.conjugate(U1_init)
            V1_init = V1_light[-1] + 0.01 * (np.random.rand(Nx) + 1j * np.random.rand(Nx))
            V2_init = np.conjugate(V1_init)