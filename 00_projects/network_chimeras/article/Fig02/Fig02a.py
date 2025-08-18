import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate
import os
import time
import datetime
from scipy.sparse import diags, hstack, vstack, identity
from scipy.linalg import expm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

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

def time_propagator(method, I, J, Q, dt):
    if method == "I_Jdt_RK4":
        k1 = J @ Q
        k2 = J @ (Q + 0.5 * dt * k1)
        k3 = J @ (Q + 0.5 * dt * k2)
        k4 = J @ (Q + dt * k3)
        return Q + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    elif method == "exp":
        return expm(J * dt) @ Q
    elif method == "I_Jdt":
        return (I + J * dt) @ Q

def jacobians_FD(eq, fields, t_grid, x_grid, y_grid, parameters, operators):
    if eq == 'duffing':
        U_1 = fields[0]
        alpha, mu, gamma, k, w = parameters  # Unpacking parameters
        DD = operators[0]

        N = len(U_1)

        j_11 = np.zeros(N)
        j_12 = np.ones(N)
        j_21 = -1 + 3 * alpha * U_1**2 - 5 * U_1**4
        j_22 = -mu * np.ones(N)

        # Use sparse matrices directly to avoid dense conversions
        J_11 = diags(j_11, 0, format='csr')
        J_12 = diags(j_12, 0, format='csr')
        J_21 = diags(j_21, 0, format='csr') + k * DD
        J_22 = diags(j_22, 0, format='csr')

        # Efficiently construct the Jacobian matrix using sparse stacking
        J_1 = hstack([J_11, J_12], format='csr')
        J_2 = hstack([J_21, J_22], format='csr')
        J = vstack([J_1, J_2], format='csr')

        return J

if __name__ == '__main__':

    project_name = '/network_chimeras/FIG02/erdos_renyi'
    disc = 'D:/'
    route = 'mnustes_science/simulation_data/FD'
    mean_degrees = [22]         #DESDE MEANDEGREE = 25 PARA ARRIBA HAY QUE USAR T TOTAL 15000, PARA ABAJO DE ESE VALOR SE PUEDE USAR T TOTAL DE 10000
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
            now = datetime.datetime.now()
            print('Hora de Término: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
            time_fin = time.time()
            print(str(time_fin - time_init) + ' seg')

            print("######## AVERAGES ########")
            t_init = 0.5 * tmax
            t_final = 1.0 * tmax
            i_0 = np.argmin(np.abs(t_light - t_init))
            i_f = np.argmin(np.abs(t_light - t_final))

            average_phase = integrate.simpson(phase_light_wraped[i_0:i_f], t_grid[i_0:i_f], axis=0) / (
                        t_grid[i_f] - t_grid[i_0])  ## AVERAGE WRAPPED PHASE
            omega_average = integrate.simpson(np.diff(phase_light, axis=0)[i_0:i_f] / dt, t_grid[i_0:i_f], axis=0) / (
                        t_grid[i_f] - t_grid[i_0])  ## AVERAGE FREQUENCY
            sigma2_phase = integrate.simpson((phase_light[i_0:i_f] - average_phase) ** 2, t_grid[i_0:i_f], axis=0) / (
                        t_grid[i_f] - t_grid[i_0])  ## FREQUENCY STANDARD DEVIATION
            sigma2 = integrate.simpson(((np.diff(phase_light, axis=0)[i_0:i_f] / dt) - omega_average) ** 2,
                                       t_grid[i_0:i_f], axis=0) / (
                                 t_grid[i_f] - t_grid[i_0])  ## FREQUENCY STANDARD DEVIATION
            average_module = integrate.simpson(module[i_0:i_f], t_grid[i_0:i_f], axis=0) / (t_grid[i_f] - t_grid[i_0])
            args = np.argsort(average_module)

            print("######## LYAPUNOV ########")
            T_00 = int(0.95 * Nt)
            T_01 = int(1.00 * Nt)
            I = np.eye(2 * Nx)
            lyap = []
            N_condit = int(2 * Nx)
            U_init = np.random.rand(2 * Nx, N_condit) - 0.5
            Q, R = np.linalg.qr(U_init)
            for i in range(T_00, T_01 - 1):
                U_1 = U_light[i, :]
                U_2 = V_light[i, :]
                J = jacobians_FD(eq, [U_1, U_2], t_light, x_grid, [0], parameters_np, operators)
                Q_new = time_propagator("I_Jdt_RK4", I, J, Q, dt)
                Q, R = np.linalg.qr(Q_new)
                Q = Q
                lyap.append(np.log(np.absolute(R.diagonal())) / dt)
            lyap = np.array(lyap)
            lyap_mean = []
            average_Nt_i = T_00
            for i in range(N_condit):
                lyap_mean.append(np.sum(lyap[:, i]))
            lyap_mean = np.array(lyap_mean) / len(lyap[:, 0])

            arg_lyap_pos = lyap_mean >= 0.0
            Q_pos = (np.sum(np.abs(Q[:N_nodes, arg_lyap_pos]) + np.abs(Q[N_nodes:, arg_lyap_pos]), axis=1)) / (
                len(lyap_mean[arg_lyap_pos]))
            args_lyap_modes = np.argsort(Q_pos)

            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(111, projection='3d')

            # Scatter 3D
            ax.scatter(average_phase, average_module, Q_pos, color='k', s=30)

            # Etiquetas
            ax.set_xlabel(r"$\langle \theta_i \rangle$", fontsize=14)
            ax.set_ylabel(r"$\langle r_i \rangle$", fontsize=14)
            ax.set_zlabel("Chaos Contribution", fontsize=14)

            # Estética
            ax.tick_params(axis='both', labelsize=12)
            ax.view_init(elev=20, azim=135)  # ajustá esto para rotar la vista

            plt.tight_layout()
            plt.show()
