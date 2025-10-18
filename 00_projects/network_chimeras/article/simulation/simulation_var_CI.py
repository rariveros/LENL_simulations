import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate
import os
import time
import datetime
from scipy.sparse import diags, hstack, vstack, identity
from scipy.linalg import expm
import gc

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
        k_2 = equations_FD(eq, old_fields + 0.5 * dt * k_1, t_grid[i] + 0.5 * dt, x_grid, y_grid, parameters, operators)
        k_3 = equations_FD(eq, old_fields + 0.5 * dt * k_2, t_grid[i] + 0.5 * dt, x_grid, y_grid, parameters, operators)
        k_4 = equations_FD(eq, old_fields + dt * k_3, t_grid[i] + dt, x_grid, y_grid, parameters, operators)
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

    project_name = '/network_chimeras/data/FIG03a/erdos_renyi'
    disc = 'G:/'
    route = 'My Drive/02. Académico/Investigación/Proyectos/[2025] Network Chimeras'
    Ks = [0.017] #np.arange(0.02, 0.025, 0.0005)         #DESDE MEANDEGREE = 25 PARA ARRIBA HAY QUE USAR T TOTAL 15000, PARA ABAJO DE ESE VALOR SE PUEDE USAR T TOTAL DE 10000
    samples = np.arange(34, 50)

    ########### NETWORK PARAMETERS ###########
    mean_degree = 20
    n = 501
    p = mean_degree / (n - 1)
    graph = erdos_renyi_graph(n, p)
    adj_matrix = np.loadtxt(r"G:\My Drive\02. Académico\Investigación\Proyectos\[2025] Network Chimeras\network_chimeras\data\FIG03b\erdos_renyi\k=0.0100\mean_degree=20.00\sample=0.00\Adj_matrix.txt", delimiter=',') #nx.adjacency_matrix(graph).toarray()
    degrees = np.sum(adj_matrix, axis=1)
    D = np.diag(degrees)

    # Laplaciana
    L_dense = D - adj_matrix
    laplacian_matrix = L_dense
    #laplacian_matrix = nx.laplacian_matrix(graph).tocsc()
    #L_dense = laplacian_matrix.toarray()
    for k in Ks:
        print("########### k = " + str(k) + " ###########")
        for sample in samples:
            print("########### Sample = " + str(sample) + " ###########")
            now = datetime.datetime.now()
            print('Hora de Inicio: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
            time_init = time.time()

            ########### DEFINING PARAMETERS ###########
            alpha = 0.4                                         # NONLINEAR COEFFICIENT
            mu = 0.1                                            # DISSIPATION
            gamma = 2.90  # 2.7  2.90                           # DRIVE STRENGTH
            #k = 0.015  # 0.4216 #0.028                          # COUPLING
            w = 0.7                                             # FORCING FREQUENCY
            eq = 'duffing'                                      # OSCILLATOR TYPE
            t_rate = 5                                          # SAVING TIME SIMULATION

            ########### PREPARING SIMULATION ########
            N_nodes = n
            [tmin, tmax, dt] = [0, 10000, 0.05]
            t_grid = np.arange(tmin, tmax + dt, dt)             # TEMPORAL GRID DEFINITION
            [xmin, xmax, dx] = [0, N_nodes, 1]
            x_grid = np.arange(xmin, xmax, dx)                  # SPATIAL GRID DEFINITION
            T = tmax
            Nt = t_grid.shape[0]
            Nx = x_grid.shape[0]

            U_init = 1.0 * np.ones(Nx)
            initial_quimera = np.random.randint(0, N_nodes)
            arg_chimera = [initial_quimera]                     # INITIAL QUIMERA INDEX
            for i in arg_chimera:
                U_init[i] = 2.0
            U_init = U_init + 0.0 * (np.random.rand(Nx) - 0.5)
            V_init = 0.0 * (np.random.rand(Nx) - 0.5)

            operators = [laplacian_matrix]
            fields_init = [U_init, V_init]
            grids = [t_grid, x_grid, 0]
            parameters_np = np.array([alpha, mu, gamma, k, w])

            ########### SIMULATION & DATA PACKING ########
            final_fields, fields_history, time_grid = RK4_FD(eq, fields_init, parameters_np, grids, dt, Nt, operators, t_rate)

            time_grid = np.array(time_grid)
            t_init = 9000
            t_final = 10000
            i_0 = np.argmin(np.abs(time_grid - t_init))
            i_f = np.argmin(np.abs(time_grid - t_final))

            U_light = np.array(fields_history)[i_0:i_f, 0]
            V_light = np.array(fields_history)[i_0:i_f, 1]

            del fields_history
            gc.collect()

            phase_light_wraped = np.arctan2(V_light, U_light)
            phase_light = np.unwrap(phase_light_wraped, axis=0)
            t_light = np.array(time_grid[i_0:i_f])
            t_grid = t_light
            module = np.sqrt(U_light ** 2 + V_light ** 2)

            del U_light, V_light
            gc.collect()

            i_0 = 0
            i_f = -1

            average_phase = integrate.simpson(phase_light_wraped[i_0:i_f], t_grid[i_0:i_f], axis=0) / (
                        t_grid[i_f] - t_grid[i_0])  ## AVERAGE WRAPPED PHASE
            average_module = integrate.simpson(module[i_0:i_f], t_grid[i_0:i_f], axis=0) / (t_grid[i_f] - t_grid[i_0])

            args = np.argsort(average_module)

            sorted_module = average_module[args]
            sorted_phase = average_phase[args]
            X_sorted = x_grid[args]

            file = disc + route + project_name
            k_str = f"{k:.{4}f}"
            mean_degree_str = f"{mean_degree:.{2}f}"
            sample_str = f"{sample:.{2}f}"
            subfile = "/k=" + k_str + "/mean_degree=" + mean_degree_str + "/sample=" + sample_str
            if not os.path.exists(file + subfile):
                os.makedirs(file + subfile)

            np.savetxt(file + subfile + '/module_mean.txt', sorted_module, delimiter=',')
            np.savetxt(file + subfile + '/phase_mean.txt', sorted_phase, delimiter=',')
            np.savetxt(file + subfile + '/X_sorted.txt', X_sorted, delimiter=',')
            np.savetxt(file + subfile + '/params.txt', parameters_np, delimiter=',')
            np.savetxt(file + subfile + '/Adj_matrix.txt', adj_matrix, delimiter=',')

            t_init = 9950
            t_final = 10000
            i_0 = np.argmin(np.abs(t_light - t_init))
            i_f = np.argmin(np.abs(t_light - t_final))

            fig, (ax01, ax02) = plt.subplots(2, 1, figsize=(7, 6))
            cax_01 = ax01.pcolormesh(x_grid, t_light[i_0:i_f], module[i_0:i_f, :], cmap="turbo", shading='auto')
            cbar = fig.colorbar(cax_01)
            cbar.ax.tick_params(labelsize=13)
            cbar.set_label('$r(t)$', rotation=0, size=20, labelpad=-50, y=1.1)
            ax01.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True,
                             labelright=False)
            ax01.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False,
                             labelbottom=True)
            ax01.set_xlabel("$i$", fontsize=20)
            ax01.set_ylabel("$t$", fontsize=20)

            cax_02 = ax02.pcolormesh(x_grid, t_light[i_0:i_f], module[i_0:i_f, args], cmap="turbo", shading='auto')
            cbar = fig.colorbar(cax_02)
            cbar.ax.tick_params(labelsize=13)
            cbar.set_label('$r(t)$', rotation=0, size=20, labelpad=-50, y=1.1)
            ax02.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True,
                             labelright=False)
            ax02.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False,
                             labelbottom=True)
            ax02.set_xlabel("$i$", fontsize=20)
            ax02.set_ylabel("$t$", fontsize=20)
            plt.savefig(file + subfile + '/spatiotemporal.png', dpi=200)

            now = datetime.datetime.now()
            print('Hora de Término: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
            time_fin = time.time()
            print(str(time_fin - time_init) + ' seg')