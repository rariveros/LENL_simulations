import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate
import os
import time
import datetime
from scipy.sparse import diags, hstack, vstack, identity
from scipy.linalg import expm

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

    project_name = '/network_chimeras/FIG01/erdos_renyi'
    disc = 'D:/'
    route = 'mnustes_science/simulation_data/FD'
    mean_degrees = [18, 26]         #DESDE MEANDEGREE = 25 PARA ARRIBA HAY QUE USAR T TOTAL 15000, PARA ABAJO DE ESE VALOR SE PUEDE USAR T TOTAL DE 10000
    samples = [0]
    for mean_degree in mean_degrees:
        print("########### MEAN DEGREE = " + str(mean_degree) + " ###########")
        for sample in samples:
            print("########### Sample = " + str(sample) + " ###########")
            now = datetime.datetime.now()
            print('Hora de Inicio: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
            time_init = time.time()

            ########### NETWORK PARAMETERS ###########
            n = 501
            p = mean_degree / (n - 1)
            graph = erdos_renyi_graph(n, p)
            adj_matrix = nx.adjacency_matrix(graph).toarray()
            laplacian_matrix = nx.laplacian_matrix(graph).tocsc()
            L_dense = laplacian_matrix.toarray()

            ########### DEFINING PARAMETERS ###########
            alpha = 0.4                                         # NONLINEAR COEFFICIENT
            mu = 0.1                                            # DISSIPATION
            gamma = 2.90  # 2.7  2.90                           # DRIVE STRENGTH
            k = 0.015  # 0.4216 #0.028                          # COUPLING
            w = 0.7                                             # FORCING FREQUENCY
            eq = 'duffing'                                      # OSCILLATOR TYPE
            t_rate = 1                                          # SAVING TIME SIMULATION

            ########### PREPARING SIMULATION ########
            N_nodes = n
            [tmin, tmax, dt] = [0, 5000, 0.05]
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

            file = disc + route + project_name
            k_str = f"{k:.{4}f}"
            mean_degree_str = f"{mean_degree:.{2}f}"
            sample_str = f"{sample:.{2}f}"
            subfile = "/k=" + k_str + "/mean_degree=" + mean_degree_str
            if not os.path.exists(file + subfile):
                os.makedirs(file + subfile)


            np.savetxt(file + subfile + '/U.txt', U_light[::10], delimiter=',')
            np.savetxt(file + subfile + '/V.txt', V_light[::10], delimiter=',')
            np.savetxt(file + subfile + '/T.txt', t_light[::10], delimiter=',')
            np.savetxt(file + subfile + '/X.txt', x_grid, delimiter=',')
            np.savetxt(file + subfile + '/params.txt', parameters_np, delimiter=',')
            np.savetxt(file + subfile + '/L.txt', L_dense, delimiter=',')

            N_init = int(0.9 * Nt)
            np.savetxt(file + subfile + '/U_lyap.txt', U_light[N_init:, :], delimiter=',')
            np.savetxt(file + subfile + '/V_lyap.txt', V_light[N_init:, :], delimiter=',')
            np.savetxt(file + subfile + '/T_lyap.txt', t_light[N_init:], delimiter=',')

            now = datetime.datetime.now()
            print('Hora de Término: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
            time_fin = time.time()
            print(str(time_fin - time_init) + ' seg')


