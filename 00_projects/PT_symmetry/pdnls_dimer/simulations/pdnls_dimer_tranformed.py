import matplotlib.pyplot as plt
from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':

    # Definiendo parámetros
    project_name = '/pdnlS_dimer_tranformed_test'
    disc = 'D:/'
    route = 'mnustes_science/simulation_data/FD'
    eq = 'pdnlS_dimer_tranformed'
    t_rate = 1

    nu = 0.1
    mu = 0.1
    Gammas = np.arange(0.18, 0.30, 0.005)
    K = [0.05]#np.arange(0.00, 0.51, 0.02)
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
    U_init = 0.01 * (np.random.rand(Nx) + 1j * np.random.rand(Nx))[0]
    V_init = 0.0 * (np.random.rand(Nx) + 1j * np.random.rand(Nx))[0]
    THETA = []
    PHI = []
    DIFF = []
    SUM = []
    R = []
    P = []
    GAMMAS = []
    for gamma in Gammas:
        for k in K:
            Delta_str = f"{nu:.{4}f}"
            gamma_str = f"{mu:.{4}f}"
            Omega_str = f"{gamma:.{4}f}"
            k_str = f"{k:.{4}f}"

            print("####### " + Omega_str + " #######")

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
            if gamma > 0.2:
                difference_phase = np.unwrap(np.angle(U_light))[-1]  - np.unwrap(np.angle(V_light))[-1]
                sum_phase = np.unwrap(np.angle(U_light))[-1] + np.unwrap(np.angle(V_light))[-1]
                THETA.append(np.unwrap(np.angle(U_light))[-1])
                PHI.append(np.unwrap(np.angle(V_light))[-1])
                DIFF.append(difference_phase)
                SUM.append(sum_phase)
                R.append(np.abs(U_light)[-1])
                P.append(np.abs(V_light)[-1])
                GAMMAS.append(gamma)
            print("R/P calculated:" + str((np.abs(U_light[-1]) / np.abs(V_light[-1])) ** 2))
            print("parameters calculated:" + str((gamma + mu) / (gamma - mu)))

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

            fig, (ax1, ax2) = plt.subplots(2)
            ax1.plot(t_light, np.abs(U_light), c="k")
            ax1.set_ylabel('$|u|$', size='20')
            ax1.set_xlim([ti * t_light[-1], t_light[-1]])
            ax1.grid(linestyle='--', alpha=0.5)

            ax2.plot(t_light, np.abs(V_light), c="k")
            ax2.set_xlabel('$t$', size='20')
            ax2.set_ylabel('$|v|$', size='20')
            ax2.set_xlim([ti * t_light[-1], t_light[-1]])
            ax2.grid(linestyle='--', alpha=0.5)
            plt.savefig(file + subfile + "/timeseries_mod.png", dpi=300)
            plt.close()

            fig, (ax1, ax2) = plt.subplots(2)
            ax1.plot(t_light, np.unwrap(np.angle(U_light)), c="k")
            ax1.set_ylabel('$\\textrm{arg} u$', size='20')
            ax1.set_xlim([ti * t_light[-1], t_light[-1]])
            ax1.grid(linestyle='--', alpha=0.5)

            ax2.plot(t_light, np.unwrap(np.angle(V_light)), c="k")
            ax2.set_xlabel('$t$', size='20')
            ax2.set_ylabel('$\\textrm{arg} v$', size='20')
            ax2.set_xlim([ti * t_light[-1], t_light[-1]])
            ax2.grid(linestyle='--', alpha=0.5)
            plt.savefig(file + subfile + "/timeseries_angle.png", dpi=300)
            plt.close()

            fig, (ax1, ax2) = plt.subplots(2)
            ax1.plot(t_light, np.unwrap(np.angle(U_light)) + np.unwrap(np.angle(V_light)), c="k")
            ax1.set_ylabel('$\\theta + \phi$', size='20')
            ax1.set_xlim([ti * t_light[-1], t_light[-1]])
            ax1.grid(linestyle='--', alpha=0.5)

            ax2.plot(t_light, np.unwrap(np.angle(U_light))  - np.unwrap(np.angle(V_light)), c="k")
            ax2.set_xlabel('$t$', size='20')
            ax2.set_ylabel('$\\theta - \phi$', size='20')
            ax2.set_xlim([ti * t_light[-1], t_light[-1]])
            ax2.grid(linestyle='--', alpha=0.5)
            plt.savefig(file + subfile + "/timeseries_relative_angles.png", dpi=300)
            plt.close()

            U_init = U_light[-1] #+ 0.01 * (np.random.rand(Nx) + 1j * np.random.rand(Nx))
            V_init = V_light[-1] #+ 0.01 * (np.random.rand(Nx) + 1j * np.random.rand(Nx))
    THETA_wrapped = (np.array(THETA) + np.pi) % (2 * np.pi) - np.pi
    PHI_wrapped = (np.array(PHI) + np.pi) % (2 * np.pi) - np.pi

    u_abs2 = ((Gammas + mu) / Gammas) * (-nu + np.sqrt(Gammas ** 2 - mu ** 2))
    v_abs2 = u_abs2 * (Gammas - mu) / (Gammas + mu)

    fig, ((ax0a, ax0b), (ax1, ax2), (ax3, ax4)) = plt.subplots(3, 2, figsize=(8, 8))
    ax0a.scatter(GAMMAS, THETA_wrapped, c="k")
    ax0a.set_ylabel('$\\theta$', size=18)
    ax0a.set_xlabel("$\gamma$", size=18)
    ax0a.grid(linestyle='--', alpha=0.5)

    ax0b.scatter(GAMMAS, PHI_wrapped, c="k")
    ax0b.set_ylabel('$\phi$', size=18)
    ax0b.set_xlabel("$\gamma$", size=18)
    ax0b.grid(linestyle='--', alpha=0.5)

    ax1.scatter(GAMMAS, SUM, c="k")
    ax1.set_ylabel('$\\theta + \phi$', size=18)
    ax1.set_xlabel("$\gamma$", size=18)
    ax1.grid(linestyle='--', alpha=0.5)

    ax2.scatter(GAMMAS, DIFF, c="k")
    ax2.set_ylabel('$\\theta - \phi$', size=18)
    ax2.set_xlabel("$\gamma$", size=18)
    ax2.grid(linestyle='--', alpha=0.5)

    ax3.scatter(GAMMAS, R, c="k")
    ax3.plot(Gammas, np.sqrt(u_abs2), color='b')
    ax3.set_ylabel('$R$', size=18)
    ax3.set_xlabel("$\gamma$", size=18)
    ax3.grid(linestyle='--', alpha=0.5)

    ax4.scatter(GAMMAS, P, c="k")
    ax4.plot(Gammas, np.sqrt(v_abs2), color='r')
    ax4.set_ylabel('$P$', size=18)
    ax4.set_xlabel("$\gamma$", size=18)
    ax4.grid(linestyle='--', alpha=0.5)

    plt.tight_layout()
    up_directory = file + "/nu=" + Delta_str + "/mu=" + gamma_str + "/k=" + k_str
    plt.savefig(up_directory + "/variables.png", dpi=300)
    plt.close()