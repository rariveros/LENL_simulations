import matplotlib.pyplot as plt
from functions import *
from back_process import *
from time_integrators import *

from numpy import correlate

if __name__ == '__main__':
    # Hay que barrer g y k, hacer espacio de parametros con fano factor o g.
    # Hay que samplear hartas realizaciones.
    # Hay que tener los diagramas de fase estocasticos y mean-field.
    # Definiendo parámetros

    project_name = '/coherent/phase_space'
    disc = 'D:/'
    route = 'mnustes_science/simulation_data/FD'
    eq = 'coherent_langevin'
    t_rate = 1

    Delta = 0.0
    gamma = 0.1
    Omegas = [0.2]
    K = [0.2, 0.42, 0.6] #np.arange(0.5, 0.1, -0.01) # [0.15] #
    g = 0.001

    # Definición de la grilla
    [tmin, tmax, dt] = [0, 5000, 0.1]
    [xmin, xmax, dx] = [0, 1, 1]
    t_grid = np.arange(tmin, tmax + dt, dt)
    x_grid = np.arange(xmin, xmax, dx)
    T = tmax
    Nt = t_grid.shape[0]
    Nx = x_grid.shape[0]
    operators = [0]
    ti = 0.

    sigma = 10.0
    U_init = 0.01 * (np.abs(np.random.normal(0, sigma, Nx)) + 1j * np.abs(np.random.normal(0, sigma, Nx)))
    V_init = 0.01 * (np.abs(np.random.normal(0, sigma, Nx)) + 1j * np.abs(np.random.normal(0, sigma, Nx)))

    U_init = np.concatenate([np.conjugate(U_init), U_init, np.conjugate(U_init), U_init, np.conjugate(U_init), U_init, np.conjugate(U_init), U_init, -np.conjugate(U_init), -U_init, -np.conjugate(U_init), -U_init, -np.conjugate(U_init), -U_init, -np.conjugate(U_init), -U_init])
    V_init = np.concatenate([V_init, np.conjugate(V_init), np.conjugate(V_init), V_init, -V_init, -np.conjugate(V_init), -np.conjugate(V_init), -V_init, V_init, np.conjugate(V_init), np.conjugate(V_init), V_init, -V_init, -np.conjugate(V_init), -np.conjugate(V_init), -V_init,])

    Re1 = np.real(U_init.flatten())
    Im1 = np.imag(U_init.flatten())
    Re2 = np.real(V_init.flatten())
    Im2 = np.imag(V_init.flatten())

    plt.figure(figsize=(6, 5))
    plt.scatter(Re1, Im1, c="k")
    plt.show()
    plt.close()

    for Omega in Omegas:
        for k in K:
            Delta_str = f"{Delta:.{4}f}"
            gamma_str = f"{gamma:.{4}f}"
            Omega_str = f"{Omega:.{4}f}"
            k_str = f"{k:.{4}f}"
            g_str = f"{g:.{4}f}"

            print("####### " + k_str + " #######")

            # Empaquetamiento de parametros, campos y derivadas para integración
            fields_init = [U_init, V_init]
            grids = [t_grid, x_grid, 0]

            parameters_np = np.array([Delta, gamma, Omega, k, g])

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
            subfile = "/Delta=" + Delta_str + "/gamma=" + gamma_str + "/Omega=" + Omega_str + "/g=" + g_str + "/k=" + k_str
            if not os.path.exists(file + subfile):
                os.makedirs(file + subfile)

            lightness = 1
            cut = int(0.2 * len(time_grid))
            U = U[cut::lightness]
            V = V[cut::lightness]
            time_grid = time_grid[cut::lightness]

            Re1 = np.real(U.flatten())
            Im1 = np.imag(U.flatten())
            Re2 = np.real(V.flatten())
            Im2 = np.imag(V.flatten())

            H1, xedges1, yedges1 = np.histogram2d(Re1, Im1, bins=200, range=[[-7, 7], [-24, 24]])#, range=[[-1.8, 1.8], [-2.5, 2.5]]) #
            H2, xedges2, yedges2 = np.histogram2d(Re2, Im2, bins=200, range=[[-24, 24], [-7, 7]])#, range=[[-2.5, 2.5], [-1.8, 1.8]]) #
            np.savetxt(file + subfile + "/histogram_counts1.txt", H1)
            np.savetxt(file + subfile + "/xedges1.txt", xedges1)
            np.savetxt(file + subfile + "/yedges1.txt", yedges1)
            np.savetxt(file + subfile + "/histogram_counts2.txt", H2)
            np.savetxt(file + subfile + "/xedges2.txt", xedges2)
            np.savetxt(file + subfile + "/yedges2.txt", yedges2)

            del Re1, Re2, Im1, Im2, U, V, time_grid
            """
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4, 8))
            hist_01 = ax1.hist2d(Re1, Im1, bins=200, density=True, cmap="inferno")
            ax1.set_xlabel("$\\textrm{Re }(\\alpha_1)$")
            ax1.set_ylabel("$\\textrm{Im }(\\alpha_1)$")

            hist_01 = ax2.hist2d(Re2, Im2, bins=200, density=True, cmap="inferno")
            ax2.set_xlabel("$\\textrm{Re }(\\alpha_2)$")
            ax2.set_ylabel("$\\textrm{Im }(\\alpha_2)$")
            #plt.gca().set_aspect("equal")
            plt.savefig("phase_space.png", dpi=300)
            plt.close()

            from scipy.ndimage import gaussian_filter

            H, xedges, yedges = np.histogram2d(Re1, Im1, bins=200, density=True)
            H_smooth = gaussian_filter(H, sigma=0.5)  # sigma controla el nivel de suavizado

            fig, ax = plt.subplots()
            ax.imshow(H_smooth.T, origin='lower', cmap="inferno",
                      extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                      aspect="auto")
            ax.set_xlabel("$\\textrm{Re }(\\alpha_2)$")
            ax.set_ylabel("$\\textrm{Im }(\\alpha_2)$")
            #plt.gca().set_aspect("equal")
            plt.savefig("kde.png", dpi=300)
            plt.close()

            U_init = U_light[-1] + 0.01 * (np.random.rand(16 * Nx) + 1j * np.random.rand(16 * Nx))
            V_init = V_light[-1] + 0.01 * (np.random.rand(16 * Nx) + 1j * np.random.rand(16 * Nx))
            """