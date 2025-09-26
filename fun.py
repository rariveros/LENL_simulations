import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate
import os
import time
import datetime
import gc
from scipy import integrate
from scipy.sparse import diags, hstack, vstack, identity
from scipy.linalg import expm

###########################################################
# FUNCIONES GENERALES
###########################################################

def erdos_renyi_graph(num_nodes, probability):
    graph = nx.erdos_renyi_graph(n=num_nodes, p=probability)
    largest_cc = max(nx.connected_components(graph), key=len)
    graph = graph.subgraph(largest_cc)
    return graph

def load_results_with_args_list(path):
    """
    RETORNA
        Lista con los siguientes elementos en orden fijo:
        [0] module         -> array módulo promedio (orden natural)
        [1] phase          -> array fase promedio (orden natural)
        [2] X              -> índices naturales
        [3] args           -> índices de reordenamiento (orden por módulo)
        [4] module_sorted  -> módulo promedio reordenado
        [5] phase_sorted   -> fase promedio reordenada
        [6] X_sorted       -> índices reordenados
    """
    # Archivos esperados
    module_path = os.path.join(path, 'module_mean.txt')
    phase_path  = os.path.join(path, 'phase_mean.txt')
    x_path      = os.path.join(path, 'X.txt')
    args_path   = os.path.join(path, 'args_module.txt')

    # Cargar datos
    module = np.loadtxt(module_path, delimiter=',')
    phase  = np.loadtxt(phase_path, delimiter=',')
    X      = np.loadtxt(x_path, delimiter=',')
    args   = np.loadtxt(args_path, delimiter=',', dtype=int)

    # Aplicar reordenamiento
    module_sorted = module[args]
    phase_sorted  = phase[args]
    X_sorted      = X[args]

    results = [module, phase, X, args, module_sorted, phase_sorted, X_sorted]
    return results

###########################################################
# TIME INTEGRATORS
###########################################################
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

###########################################################
# EQUATIONS AND OPERATORS
###########################################################

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