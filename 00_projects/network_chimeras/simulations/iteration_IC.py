import networkx as nx
import random
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.sparse import diags
from scipy import integrate
from sklearn.cluster import KMeans
import h5py

# Generate a random graph using the Erdős-Rényi model
def erdos_renyi_graph(num_nodes, probability):
    graph = nx.erdos_renyi_graph(n=num_nodes, p=probability)
    # keep the largest connected component
    largest_cc = max(nx.connected_components(graph), key=len)
    graph = graph.subgraph(largest_cc)
    return graph


# MODEL SOLVE
def RK4_FD(eq, fields, parameters, grids, dt, Nt, operators, t_rate):  # implementa rouge-kutta
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


def equations_FD(eq, field_slices, t_i, x_grid, y_grid, parameters, operators):  # ecuaciones
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


def Der(D, f):  # función de diferenciación
    d_f = D @ f
    return d_f


def time_propagator(type, I, J, Q, dt):
    if type == "I_Jdt_RK4":
        k_1 = J @ Q  # np.exp(J * dt) @ Q #
        k_2 = J @ (Q + 0.5 * dt * k_1)  # np.exp(J * dt) @ (Q + 0.5 * dt * k_1) #
        k_3 = J @ (Q + 0.5 * dt * k_2)  # np.exp(J * dt) @ (Q + 0.5 * dt * k_2) #
        k_4 = J @ (Q + dt * k_3)  # np.exp(J * dt) @ (Q + dt * k_3) #
        U_new = Q + dt * (k_1 + 2 * k_2 + 2 * k_3 + k_4) / 6
    elif type == "exp":
        U_new = np.exp(J * dt) @ Q
    elif type == "I_Jdt":
        U_new = (I + J * dt) @ Q
    return U_new


def jacobians_FD(eq, fields, t_grid, x_grid, y_grid, parameters, operators):
    if eq == 'duffing':
        U_1 = fields[0]
        U_2 = fields[1]

        alpha = parameters[0]  # NONLINEAR COEFFICIENT
        mu = parameters[1]  # DISSIPATION
        gamma = parameters[2]  # 2.7                                    # DRIVE STRENGTH
        k = parameters[3]  # 0.4216                                      # COUPLING (0.42)
        w = parameters[4]

        DD = operators[0]

        j_11 = [np.zeros(len(U_1))]
        j_12 = [np.ones(len(U_1))]
        j_21 = [- 1 + 3 * alpha * U_1 ** 2 - 5 * U_1 ** 4]
        j_22 = [-mu * np.ones(len(U_1))]

        J_11 = diags(j_11, [0])
        J_11 = J_11.toarray()

        J_12 = diags(j_12, [0])
        J_12 = J_12.toarray()

        J_21 = diags(j_21, [0]) + k * DD
        J_21 = J_21.toarray()

        J_22 = diags(j_22, [0])
        J_22 = J_22.toarray()

        J_1 = np.concatenate((J_11, J_12), axis=1)
        J_2 = np.concatenate((J_21, J_22), axis=1)
        J = np.concatenate((J_1, J_2), axis=0)
    return J


def classification_by_moldule_threshold(average_module, power_threshold):
    arg_chimeras = average_module < 1.5
    arg_sync = average_module >= 1.5

    return arg_chimeras, arg_sync


def compute_lyapunov_exponents(Nt, Nx, U_light, V_light, eq, t_grid, x_grid, parameters_np, operators, dt):
    T_00 = int(0.9 * Nt)
    T_01 = int(1.00 * Nt)
    I = np.eye(2 * Nx)
    lyap_01 = []
    N_condit = int(2 * Nx)
    U_init = np.random.rand(2 * Nx, N_condit) - 0.5
    Q, R = np.linalg.qr(U_init)

    for i in range(T_00, T_01 - 1):
        U_1 = U_light[i, :]
        U_2 = V_light[i, :]
        J = jacobians_FD(eq, [U_1, U_2], t_grid, x_grid, [0], parameters_np, operators)
        Q_new = time_propagator("I_Jdt_RK4", I, J, Q, dt)
        Q, R = np.linalg.qr(Q_new)
        lyap_01.append(np.log(np.absolute(R.diagonal())) / dt)

    lyap_01 = np.array(lyap_01)
    lyap_mean_01 = np.array([np.sum(lyap_01[:, i]) for i in range(N_condit)]) / len(lyap_01[:, 0])

    return lyap_mean_01, Q, N_condit


def lyap_contribution(lyap_mean_01, Q, N_nodes):
    arg_lyap_pos = lyap_mean_01 >= 0.0
    Q_pos = (np.sum(np.abs(Q[:N_nodes, arg_lyap_pos]) + np.abs(Q[N_nodes:, arg_lyap_pos]), axis=1)) / (
        len(lyap_mean_01[arg_lyap_pos]))

    return Q_pos


def classification_by_kmeans(Q_pos, average_phase):
    data = list(zip(Q_pos, average_phase))

    kmeans = KMeans(n_clusters=2)
    kmeans.fit(data)
    centroids = kmeans.cluster_centers_
    # clasification Syncronized & Chimera

    if centroids[0, 0] > centroids[1, 0]:
        arg_chimeras = kmeans.labels_ == 0
        arg_sync = kmeans.labels_ == 1
    elif centroids[1, 0] > centroids[0, 0]:
        arg_chimeras = kmeans.labels_ == 1
        arg_sync = kmeans.labels_ == 0
    else:
        raise ValueError("No classification is possible")

    return arg_chimeras, arg_sync


# Model - Parameter
if __name__ == '__main__':
    alpha = 0.4  # NONLINEAR COEFFICIENT
    mu = 0.1  # DISSIPATION
    gamma = 2.90  # 2.7  2.90                                  # DRIVE STRENGTH  # 0.4216 #0.028                                  # COUPLING (0.42)
    w = 0.7
    eq = 'duffing'
    t_rate = 1

    # Grid definition
    N_nodes = 100
    [tmin, tmax, dt] = [0, 5000, 0.05]
    t_grid = np.arange(tmin, tmax + dt, dt)  # TEMPORAL GRID DEFINITION
    [xmin, xmax, dx] = [0, N_nodes, 1]
    x_grid = np.arange(xmin, xmax, dx)  # SPATIAL GRID DEFINITION

    T = tmax
    Nt = t_grid.shape[0]
    Nx = x_grid.shape[0]

    n = 10
    p = n / N_nodes
    P = [p]
    IC = [40]
    Ks = np.arange(0.015, 0.07, 0.0025)

    AVERAGE_PHASE = {}
    OMEGA_AVERAGE = {}
    SIGMA2 = {}
    AVERAGE_MODULE = {}
    AV2_MODULE = {}
    FRAC_CHIMERA = {}
    NODE_CENTRALITY = {}

    graph = erdos_renyi_graph(N_nodes, p)
    centrality = nx.degree_centrality(graph)
    node_centrality = list(centrality.values())
    adj_matrix = nx.adjacency_matrix(graph).toarray()
    laplacian_matrix = nx.laplacian_matrix(graph).tocsc()

    for i in Ks:
        print("Pinched Node: " + str(i))# loop for every node
        # Initial Conditions
        U_init = 1.0 * np.ones(Nx)
        initial_quimera = IC[0]
        arg_chimera = [initial_quimera]  # INITIAL QUIMERA INDEX
        for j in arg_chimera:
            U_init[j] = 2.0
        U_init = U_init + 0.0 * (np.random.rand(Nx) - 0.5)
        V_init = 0.0 * np.random.rand(Nx)
        k = i
        # Empaquetamiento de parametros, campos y derivadas para integración

        operators = [laplacian_matrix]
        fields_init = [U_init, V_init]
        grids = [t_grid, x_grid, 0]
        parameters_np = np.array([alpha, mu, gamma, k, w])

        # numerical simulation
        final_fields, fields_history, time_grid = RK4_FD(eq, fields_init, parameters_np, grids, dt, Nt, operators, t_rate)

        # Reobteniendo campos y variables
        U = np.array(fields_history)[:, 0]
        V = np.array(fields_history)[:, 1]
        lightness = 1
        U_light = U[0::lightness]
        V_light = V[0::lightness]
        phase_light_wraped = np.arctan2(V_light, U_light)
        phase_light = np.unwrap(phase_light_wraped, axis=0)[0::lightness]
        t_light = np.array(time_grid[0::lightness])
        module = np.sqrt(U_light ** 2 + V_light ** 2)  ## AMPLITUDE / MODULE

        # Save results in HDF5 format
        #    np.savetxt(os.path.join(output_dir+"/sim/U_light/", f"U_light_{i}.txt"), U_light, fmt="%.6f")
        #    np.savetxt(os.path.join(output_dir+"/sim/V_light/", f"V_light_{i}.txt"), V_light, fmt="%.6f")

        ############## AVERAGES AND CHARACTERIZATION OF THE STATIONARY FINAL STATE N = (4000, 5000)##############
        t_init = 4000
        t_final = 5000

        i_0 = np.argmin(np.abs(t_light - t_init))
        i_f = np.argmin(np.abs(t_light - t_final))

        average_phase = integrate.simpson(phase_light_wraped[i_0:i_f], x=t_grid[i_0:i_f], axis=0) / (
                    t_grid[i_f] - t_grid[i_0])  ## AVERAGE WRAPPED PHASE
        omega_average = integrate.simpson(np.diff(phase_light, axis=0)[i_0:i_f] / dt, x=t_grid[i_0:i_f], axis=0) / (
                    t_grid[i_f] - t_grid[i_0])  ## AVERAGE FREQUENCY
        sigma2 = integrate.simpson(((np.diff(phase_light, axis=0)[i_0:i_f] / dt) - omega_average) ** 2, x=t_grid[i_0:i_f],
                                   axis=0) / (t_grid[i_f] - t_grid[i_0])  ## FREQUENCY STANDARD DEVIATION
        average_module = integrate.simpson(module[i_0:i_f], x=t_grid[i_0:i_f], axis=0) / (t_grid[i_f] - t_grid[i_0])
        av2_module = np.mean(average_module)

        arg_chimeras, arg_sync = classification_by_moldule_threshold(average_module, 1.5)

        frac_chimera = len(average_module[arg_chimeras]) / N_nodes
        plt.scatter(i, frac_chimera, c="k")

        AVERAGE_PHASE[i] = average_phase
        OMEGA_AVERAGE[i] = omega_average
        SIGMA2[i] = sigma2
        AVERAGE_MODULE[i] = average_module
        AV2_MODULE[i] = av2_module
        FRAC_CHIMERA[i] = frac_chimera
    plt.show()