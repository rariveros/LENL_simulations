from back_process import *
from functions import *
from hamiltonians import *


def RK4_FD(eq, fields, parameters, grids, dt, Nt, operators, t_rate):
    t_grid = grids[0]
    x_grid = grids[1]
    y_grid = grids[2]
    fields_history = []
    time_grid = []
    for i in range(Nt - 1):
        old_fields = fields
        k_1 = equations_FD(eq, old_fields, t_grid, x_grid, y_grid, parameters, operators)
        k_2 = equations_FD(eq, old_fields + 0.5 * dt * k_1, t_grid, x_grid, y_grid, parameters, operators)
        k_3 = equations_FD(eq, old_fields + 0.5 * dt * k_2, t_grid, x_grid, y_grid, parameters, operators)
        k_4 = equations_FD(eq, old_fields + dt * k_3, t_grid, x_grid, y_grid, parameters, operators)
        new_fields = old_fields + dt * (k_1 + 2 * k_2 + 2 * k_3 + k_4) / 6
        fields = new_fields
        if i % t_rate == 0:
            fields_history.append(fields)
            time_grid.append(t_grid[i])
    return fields, fields_history, time_grid


def RK4_FD_ham(eq, fields, parameters, grids, dt, Nt, operators, t_rate):
    t_grid = grids[0]
    x_grid = grids[1]
    y_grid = grids[2]
    fields_history = []
    time_grid = []
    for i in range(Nt - 1):
        old_fields = fields
        k_1 = equations_FD(eq, old_fields, t_grid, x_grid, y_grid, parameters, operators)
        k_2 = equations_FD(eq, old_fields + 0.5 * dt * k_1, t_grid, x_grid, y_grid, parameters, operators)
        k_3 = equations_FD(eq, old_fields + 0.5 * dt * k_2, t_grid, x_grid, y_grid, parameters, operators)
        k_4 = equations_FD(eq, old_fields + dt * k_3, t_grid, x_grid, y_grid, parameters, operators)
        new_fields = old_fields + dt * (k_1 + 2 * k_2 + 2 * k_3 + k_4) / 6
        fields = new_fields
        hamiltonian = hamiltonians_FD(eq, fields, t_grid, x_grid, y_grid, parameters, operators)
        if i % t_rate == 0:
            fields_history.append(hamiltonian)
            time_grid.append(t_grid[i])
    return fields, fields_history, time_grid


def RK4_ode(eq, vector, t_grid, parameters, dt, Nt):
    vector_history = []
    for i in range(Nt - 1):
        old_vector = vector
        k_1 = equations_ode(eq, old_vector, t_grid[i], parameters)
        k_2 = equations_ode(eq, old_vector + 0.5 * dt * k_1, t_grid[i], parameters)
        k_3 = equations_ode(eq, old_vector + 0.5 * dt * k_2, t_grid[i], parameters)
        k_4 = equations_ode(eq, old_vector + dt * k_3, t_grid[i], parameters)
        new_vector = old_vector + dt * (k_1 + 2 * k_2 + 2 * k_3 + k_4) / 6
        vector = new_vector
        if i % 10 == 0:
            vector_history.append(vector)
    vector_history = np.array(vector_history)
    return vector, vector_history