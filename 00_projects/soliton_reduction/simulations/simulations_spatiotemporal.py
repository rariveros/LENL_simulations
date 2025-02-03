import matplotlib.pyplot as plt
from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':

    # Definiendo parámetros
    project_name = '/soliton_reduced'
    disc = 'D:/'
    route = 'mnustes_science/simulation_data/FD'
    eq = 'reduced_soliton'
    t_rate = 1

    [alpha, beta, mu, sigma, gamma] = [6.524, 1, 0.075, 15, 0.18]
    beta_adim = 0.004811649356064012
    alpha_str = f"{alpha:.{3}f}"
    beta_str = f"{beta:.{3}f}"
    mu_str = f"{mu:.{3}f}"
    sigma_str = f"{sigma:.{2}f}"
    gamma_str = f"{gamma:.{3}f}"

    working_directory = 'D:/mnustes_science/simulation_data/FD/soliton_reduced/alpha=' + alpha_str + '/beta=' + beta_str + '/mu=' + mu_str
    directories = ["nu=-0.050", "nu=-0.150"]
    Zs = []
    ics = [-4, 4]
    for directory in directories:
        for ic in ics:
            dir = working_directory + "/" + directory
            D1 = np.loadtxt(dir + '/gamma=' + gamma_str + '/sigma=' + sigma_str + '/D1.txt', delimiter=',')
            D2 = np.loadtxt(dir + '/gamma=' + gamma_str + '/sigma=' + sigma_str + '/D2.txt', delimiter=',')
            parameters = np.loadtxt(dir + '/gamma=' + gamma_str + '/sigma=' + sigma_str + '/parameters.txt', delimiter=',')

            ti = 0.

            # Definición de la grilla
            [tmin, tmax, dt] = [0, 300, 1]
            [xmin, xmax, dx] = [0, 1, 1]
            t_grid = np.arange(tmin, tmax + dt, dt)
            x_grid = np.arange(xmin, xmax, dx)
            T = tmax
            Nt = t_grid.shape[0]
            Nx = x_grid.shape[0]

            # Initial Conditions
            U_init = ic * np.ones(Nx)
            V_init = ic * np.ones(Nx)
            operators = [D1, D2]

            # Empaquetamiento de parametros, campos y derivadas para integración
            fields_init = [U_init, V_init]
            grids = [t_grid, x_grid, 0]
            parameters_np = parameters #[alpha, beta, mu, nu, sigma, gamma]
            [alpha, beta, mu, nu, sigma, gamma] = parameters
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
            subfile = "/test"
            if not os.path.exists(file + subfile):
                os.makedirs(file + subfile)

            np.savetxt(file + subfile + '/U.txt', U, delimiter=',')
            np.savetxt(file + subfile + '/V.txt', V, delimiter=',')
            np.savetxt(file + subfile + '/parameters.txt', parameters_np, delimiter=',')
            np.savetxt(file + subfile + '/T.txt', t_grid, delimiter=',')

            lightness = 1
            U_light = U[0::lightness]
            V_light = V[0::lightness]
            t_light = time_grid[0::lightness]

            dx = 0.2
            x_grid = np.arange(-30, 30, dx)
            PHI_01 = []
            PHI_02 = []
            for i in range(len(t_light)):
                [phi_01, phi_02, Dphi_01, Dphi_02] = phis(alpha, beta, nu, mu, gamma, sigma, U_light[i], V_light[i], x_grid, dx)
                PHI_01.append(phi_01)
                PHI_02.append(phi_02)
            PHI_01 = np.array(PHI_01)
            PHI_02 = np.array(PHI_02)

            Z = np.abs(PHI_01 + 1j * PHI_02) / np.sqrt(beta_adim)
            Zs.append(Z)
    ls_a = 18
    ls_b = 15
    ls_c = 20

    # Determine global vmin and vmax for all plots
    vmin = 0
    vmax = max(np.max(np.real(Zs[0])), np.max(np.real(Zs[1])), np.max(np.real(Zs[2])), np.max(np.real(Zs[3])))

    fig, axs = plt.subplots(2, 2, figsize=(6, 5))  # Create a 2x2 grid of subplots
    ax01, ax02, ax03, ax04 = axs[0, 0], axs[0, 1], axs[1, 0], axs[1, 1]

    # Plot data on each subplot with consistent vmin and vmax
    cax_01 = ax01.pcolormesh(x_grid, t_light, np.real(Zs[0]), cmap=parula_map, shading='auto', vmin=vmin, vmax=vmax)
    ax01.tick_params(axis="y", direction="in", labelsize=ls_a, left=True, right=True, labelleft=True, labelright=False)
    ax01.tick_params(axis="x", direction="in", labelsize=ls_a, top=True, bottom=True, labeltop=False, labelbottom=False)
    ax01.set_ylabel("$t/T$", fontsize=ls_c)

    cax_02 = ax02.pcolormesh(x_grid, t_light, np.real(Zs[1]), cmap=parula_map, shading='auto', vmin=vmin, vmax=vmax)
    ax02.tick_params(axis="y", direction="in", labelsize=ls_a, left=True, right=True, labelleft=False, labelright=False)
    ax02.tick_params(axis="x", direction="in", labelsize=ls_a, top=True, bottom=True, labeltop=False, labelbottom=False)

    cax_03 = ax03.pcolormesh(x_grid, t_light, np.real(Zs[2]), cmap=parula_map, shading='auto', vmin=vmin, vmax=vmax)
    ax03.tick_params(axis="y", direction="in", labelsize=ls_a, left=True, right=True, labelleft=True, labelright=False)
    ax03.tick_params(axis="x", direction="in", labelsize=ls_a, top=True, bottom=True, labeltop=False, labelbottom=True)
    ax03.set_xlabel("$x\ \\textrm{(mm)}$", fontsize=ls_c)
    ax03.set_ylabel("$t/T$", fontsize=ls_c)

    cax_04 = ax04.pcolormesh(x_grid, t_light, np.real(Zs[3]), cmap=parula_map, shading='auto', vmin=vmin, vmax=vmax)
    ax04.tick_params(axis="y", direction="in", labelsize=ls_a, left=True, right=True, labelleft=False, labelright=False)
    ax04.tick_params(axis="x", direction="in", labelsize=ls_a, top=True, bottom=True, labeltop=False, labelbottom=True)
    ax04.set_xlabel("$\\varphi(x, \\xi_1)$", fontsize=ls_c)

    # Add a single colorbar connected to all plots
    cbar = fig.colorbar(cax_03, ax=axs, orientation='horizontal', fraction=0.05, pad=0.05, aspect=40, shrink=1.07, location='top')
    cbar.ax.tick_params(labelsize=ls_b, direction='in', top=True, bottom=False, labeltop=True, labelbottom=False)

    # Adjust layout to ensure subplots are aligned with the colorbar
    plt.subplots_adjust(wspace=0.3, hspace=0.0)  # Adjust spacing between plots
    plt.tight_layout(rect=[0, 0, 1, 0.85])  # Reserve space for the colorbar above the plots

    # Save the final plot
    plt.savefig(file + subfile + '/realfields_aligned_colorbar.png', dpi=300)
    plt.close()