import networkx as nx
import random
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import diags
from scipy.sparse import diags, hstack, vstack, identity
from scipy.linalg import expm
from back_process import filtro_array
from scipy import integrate
import os
import time
import datetime

def erdos_renyi_graph(num_nodes, probability):
    graph = nx.erdos_renyi_graph(n=num_nodes, p=probability)
    # keep the largest connected component
    largest_cc = max(nx.connected_components(graph), key=len)
    graph = graph.subgraph(largest_cc)
    return graph

def RK4_FD(eq, fields, parameters, grids, dt, Nt, operators, t_rate): #implementa rouge-kutta
    t_grid = grids[0]
    x_grid = grids[1]
    y_grid = grids[2]
    fields_history = []
    time_grid = []
    for i in range(Nt - 1):
        old_fields = fields
        k_1 = equations_FD(eq, old_fields, t_grid[i], x_grid, y_grid, parameters, operators)
        k_2 = equations_FD(eq, old_fields + 0.5 * dt * k_1, t_grid[i], x_grid, y_grid, parameters, operators)
        k_3 = equations_FD(eq, old_fields + 0.5 * dt * k_2, t_grid[i], x_grid, y_grid, parameters, operators)
        k_4 = equations_FD(eq, old_fields + dt * k_3, t_grid[i], x_grid, y_grid, parameters, operators)
        new_fields = old_fields + dt * (k_1 + 2 * k_2 + 2 * k_3 + k_4) / 6
        fields = new_fields
        if i % t_rate == 0:
            fields_history.append(fields)
            time_grid.append(t_grid[i])
    return fields, fields_history, time_grid

def equations_FD(eq, field_slices, t_i, x_grid, y_grid, parameters, operators): #ecuaciones
    if eq == 'duffing':
        U = field_slices[0]
        V = field_slices[1]

        alpha = parameters[0]
        mu = parameters[1]
        gamma = parameters[2]
        k = parameters[3]
        w = parameters[4]
        DD = operators[0]

        ddU = DD @ U

        F = V
        G = - U + alpha * U ** 3 - U ** 5 - mu * V + gamma * np.cos(w * t_i) + k * ddU

        fields = np.array([F, G])
    return fields

def Der(D, f): #función de diferenciación
    d_f = D @ f
    return d_f

if __name__ == '__main__':

    project_name = '/network_chimeras/erdos_renyi'
    disc = 'D:/'
    route = 'mnustes_science/simulation_data/FD'
    mean_degrees = np.arange(25, 33, 1)
    samples = [0, 1, 2]
    for mean_degree in mean_degrees:
        print("########### MEAN DEGREE = " + str(mean_degree) + " ###########")
        for sample in samples:
            print("########### Sample = " + str(sample) + " ###########")
            now = datetime.datetime.now()
            print('Hora de Inicio: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
            time_init = time.time()

            ########### NETWORK PARAMETERS ###########
            n = 501
            #mean_degree = 30
            p = mean_degree / (n - 1)
            graph = erdos_renyi_graph(n, p)
            adj_matrix = nx.adjacency_matrix(graph).toarray()
            laplacian_matrix = nx.laplacian_matrix(graph).tocsc()

            ########### DEFINING PARAMETERS ###########
            alpha = 0.4                                         # NONLINEAR COEFFICIENT
            mu = 0.1                                            # DISSIPATION
            gamma = 2.90  # 2.7  2.90                           # DRIVE STRENGTH
            k = 0.014  # 0.4216 #0.028                          # COUPLING
            w = 0.7                                             # FORCING FREQUENCY
            eq = 'duffing'                                      # OSCILLATOR TYPE
            t_rate = 1                                          # SAVING TIME SIMULATION

            ########### PREPARING SIMULATION ########
            N_nodes = n
            [tmin, tmax, dt] = [0, 15000, 0.05]
            t_grid = np.arange(tmin, tmax + dt, dt)             # TEMPORAL GRID DEFINITION
            [xmin, xmax, dx] = [0, N_nodes, 1]
            x_grid = np.arange(xmin, xmax, dx)                  # SPATIAL GRID DEFINITION
            T = tmax
            Nt = t_grid.shape[0]
            Nx = x_grid.shape[0]

            U_init = 1.0 * np.ones(Nx)
            initial_quimera = 18
            arg_chimera = [initial_quimera]                     # INITIAL QUIMERA INDEX
            for i in arg_chimera:
                U_init[i] = 2.0
            U_init = U_init + 0.0 * (np.random.rand(Nx) - 0.5)
            V_init = 0.0 * np.random.rand(Nx)

            operators = [laplacian_matrix]
            fields_init = [U_init, V_init]
            grids = [t_grid, x_grid, 0]
            parameters_np = np.array([alpha, mu, gamma, k, w])

            ########### SIMULATION & DATA PACKING ########
            final_fields, fields_history, time_grid = RK4_FD(eq, fields_init, parameters_np, grids, dt, Nt, operators, t_rate)

            U = np.array(fields_history)[:, 0]
            V = np.array(fields_history)[:, 1]
            phase = np.arctan2(V, U)
            lightness = 1
            U_light = U[0::lightness]
            V_light = V[0::lightness]
            phase_light_wraped = phase[0::lightness]
            phase_light = np.unwrap(phase_light_wraped, axis=0)[0::lightness]
            t_light = np.array(time_grid[0::lightness])
            module = np.sqrt(U_light ** 2 + V_light ** 2)

            ############### SAVING PLOTS & DATA ###############
            t_init = 14900
            t_final = 15000
            i_0 = np.argmin(np.abs(t_light - t_init))
            i_f = np.argmin(np.abs(t_light - t_final))
            power_threshold = 1.5

            #### SPATIOTEMPORAL DIAGRAMS ####

            file = disc + route + project_name
            k_str = f"{k:.{4}f}"
            mean_degree_str = f"{mean_degree:.{2}f}"
            sample_str = f"{sample:.{2}f}"
            subfile = "/k=" + k_str + "/mean_degree=" + mean_degree_str + "/sample=" + sample_str
            if not os.path.exists(file + subfile):
                os.makedirs(file + subfile)

            fig, (ax01, ax02) = plt.subplots(2, 1, figsize=(5, 4))
            cax_01 = ax01.pcolormesh(x_grid, t_light[i_0:i_f], U_light[i_0:i_f, :], cmap="turbo", shading='auto')
            cbar = fig.colorbar(cax_01)
            cbar.ax.tick_params(labelsize=13)
            cbar.set_label('$u(t)$', rotation=0, size=20, labelpad=-50, y=1.1)
            ax01.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
            ax01.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=True)
            ax01.set_xlabel("$i$", fontsize=20)
            ax01.set_ylabel("$t$", fontsize=20)

            cax_02 = ax02.pcolormesh(x_grid, t_light[i_0:i_f], module[i_0:i_f, :], cmap="turbo", shading='auto')
            cbar = fig.colorbar(cax_02)
            cbar.ax.tick_params(labelsize=13)
            cbar.set_label('$r(t)$', rotation=0, size=20, labelpad=-50, y=1.1)
            ax02.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
            ax02.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=True)
            ax02.set_xlabel("$i$", fontsize=20)
            ax02.set_ylabel("$t$", fontsize=20)

            plt.savefig(file + subfile + '/spatiotemporal.png', dpi=200)
            plt.close()

            #### AVERAGES AND CHARACTERIZATION OF THE STATIONARY FINAL STATE ####

            t_init =14000
            t_final = 15000
            i_0 = np.argmin(np.abs(t_light - t_init))
            i_f = np.argmin(np.abs(t_light - t_final))

            average_phase = integrate.simpson(phase_light_wraped[i_0:i_f], t_grid[i_0:i_f], axis=0) / (t_grid[i_f] - t_grid[i_0])                                   ## AVERAGE WRAPPED PHASE
            omega_average = integrate.simpson(np.diff(phase_light, axis=0)[i_0:i_f] / dt, t_grid[i_0:i_f], axis=0) / (t_grid[i_f] - t_grid[i_0])                    ## AVERAGE FREQUENCY
            sigma2_phase = integrate.simpson((phase_light[i_0:i_f] - average_phase) ** 2, t_grid[i_0:i_f], axis=0) / (t_grid[i_f] - t_grid[i_0])                    ## FREQUENCY STANDARD DEVIATION
            sigma2 = integrate.simpson(((np.diff(phase_light, axis=0)[i_0:i_f] / dt) - omega_average) ** 2, t_grid[i_0:i_f], axis=0) / (t_grid[i_f] - t_grid[i_0])  ## FREQUENCY STANDARD DEVIATION
            average_module = integrate.simpson(module[i_0:i_f], t_grid[i_0:i_f], axis=0) / (t_grid[i_f] - t_grid[i_0])

            args = np.argsort(average_module)

            fig, (ax01, ax02, ax03, ax04, ax05) = plt.subplots(5, 1, figsize=(10, 10))

            ax01.scatter(x_grid, average_phase[args], color="k")
            ax01.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
            ax01.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=False)
            ax01.set_ylabel("$\\langle \\theta_i \\rangle$", fontsize=20)

            ax02.scatter(x_grid, average_module[args], color="k")
            ax02.hlines(power_threshold, x_grid[0], x_grid[-1], colors="r")
            ax02.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
            ax02.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=False)
            ax02.set_ylabel("$\\langle  r_i \\rangle$", fontsize=20)

            ax03.scatter(x_grid, omega_average[args], color="k")
            ax03.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
            ax03.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=False)
            ax03.set_ylabel("$\\langle \\omega_i \\rangle$", fontsize=20)

            ax04.scatter(x_grid, np.sqrt(sigma2[args]), color="k")
            ax04.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
            ax04.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=False)
            ax04.set_ylabel("$\\sigma_i$", fontsize=20)

            ax05.scatter(x_grid, phase_light_wraped[-1, args], color="k")
            ax05.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
            ax05.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=True)
            ax05.set_ylabel("$\\theta_i(t_f)$", fontsize=20)
            ax05.set_xlabel("$i$", fontsize=20)

            plt.savefig(file + subfile + '/characterization.png', dpi=200)
            plt.close()

            #### SPATIOTEMPORAL DYNAMICS OF QUIMERAS SORTED BY PHASE/AMPLITUDE ####
            t_init = 14900
            t_final = 15000
            i_0 = np.argmin(np.abs(t_light - t_init))
            i_f = np.argmin(np.abs(t_light - t_final))

            list_spl = np.array([nx.shortest_path_length(graph, source=18, target=i) for i in graph.nodes()])
            arg_spl = np.argsort(list_spl)

            fig, (ax01, ax02) = plt.subplots(2, 1, figsize=(5, 4))
            cax_01 = ax01.pcolormesh(x_grid, t_light[i_0:i_f], U_light[i_0:i_f, args], cmap="turbo", shading='auto')
            cbar = fig.colorbar(cax_01)
            cbar.ax.tick_params(labelsize=13)
            cbar.set_label('$u(t)$', rotation=0, size=20, labelpad=-50, y=1.1)
            ax01.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
            ax01.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=True)
            ax01.set_xlabel("$i$", fontsize=20)
            ax01.set_ylabel("$t$", fontsize=20)

            cax_02 = ax02.pcolormesh(x_grid, t_light[i_0:i_f], module[i_0:i_f, args], cmap="turbo", shading='auto')
            cbar = fig.colorbar(cax_02)
            cbar.ax.tick_params(labelsize=13)
            cbar.set_label('$r(t)$', rotation=0, size=20, labelpad=-50, y=1.1)
            ax02.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
            ax02.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=True)
            ax02.set_xlabel("$i$", fontsize=20)
            ax02.set_ylabel("$t$", fontsize=20)

            plt.savefig(file + subfile + '/spatiotemporal_sorted.png', dpi=200)
            plt.close()

            ############### SAVING PLOTS & DATA ###############
            arg_chimeras = average_module < power_threshold
            arg_sync = average_module >= power_threshold

            frac_chimera = len(average_module[arg_chimeras]) / N_nodes
            output_parameters = [frac_chimera, mean_degree]
            print("[N, k] = " + str(output_parameters))
            np.savetxt(file + subfile + '/module_mean.txt', average_module[args], delimiter=',')
            np.savetxt(file + subfile + '/phase_mean.txt', average_phase[args], delimiter=',')
            np.savetxt(file + subfile + '/X_sorted.txt', x_grid[args], delimiter=',')
            np.savetxt(file + subfile + '/parameters.txt', parameters_np, delimiter=',')
            np.savetxt(file + subfile + '/output.txt', output_parameters, delimiter=',')
            np.savetxt(file + subfile + '/T.txt', t_grid, delimiter=',')

            now = datetime.datetime.now()
            print('Hora de Término: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
            time_fin = time.time()
            print(str(time_fin - time_init) + ' seg')

