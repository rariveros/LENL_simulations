from back_process import *
from functions import *
from hamiltonians import *
from jacobians import *

"""
def RK4_FD(eq, fields, parameters, grids, dt, Nt, operators, t_rate):
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
"""

def RK4_FD(eq, fields, parameters, grids, dt, Nt, operators, t_rate):
    t_grid = grids[0]  # Time grid
    x_grid = grids[1]  # Spatial grid (x)
    y_grid = grids[2]  # Spatial grid (y)

    fields_history = []
    time_grid = []

    for i in range(Nt - 1):
        t = t_grid[i]  # Current time
        old_fields = fields.copy()  # Ensure copy to avoid overwriting

        # RK4 steps with correct time updates
        k_1 = equations_FD(eq, old_fields, t, x_grid, y_grid, parameters, operators)
        k_2 = equations_FD(eq, old_fields + 0.5 * dt * k_1, t + 0.5 * dt, x_grid, y_grid, parameters, operators)
        k_3 = equations_FD(eq, old_fields + 0.5 * dt * k_2, t + 0.5 * dt, x_grid, y_grid, parameters, operators)
        k_4 = equations_FD(eq, old_fields + dt * k_3, t + dt, x_grid, y_grid, parameters, operators)

        # Update fields using weighted sum
        new_fields = old_fields + dt * (k_1 + 2 * k_2 + 2 * k_3 + k_4) / 6
        fields = new_fields  # Update for next iteration

        # Save history every `t_rate` steps
        if i % t_rate == 0:
            fields_history.append(fields.copy())  # Store a copy to avoid referencing issues
            time_grid.append(t)

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

def RK5_FD(eq, fields, parameters, grids, dt, Nt, operators, t_rate):
    t_grid = grids[0]       # Time grid
    x_grid = grids[1]       # Spatial grid (x)
    y_grid = grids[2]       # Spatial grid (y)

    fields_history = []
    time_grid = []

    # Coefficients for Cash-Karp RK5 (Fehlberg tableau)
    a = [0, 1/5, 3/10, 3/5, 1, 7/8]
    b = [
        [],
        [1/5],
        [3/40, 9/40],
        [3/10, -9/10, 6/5],
        [-11/54, 5/2, -70/27, 35/27],
        [1631/55296, 175/512, 575/13824, 44275/110592, 253/4096]
    ]
    c = [37/378, 0, 250/621, 125/594, 0, 512/1771]  # 5th-order solution

    for i in range(Nt - 1):
        t = t_grid[i]
        old_fields = fields.copy()

        # Preallocate k1 to k6
        k = [None] * 6

        # Compute stages
        k[0] = equations_FD(eq, old_fields, t, x_grid, y_grid, parameters, operators)
        for j in range(1, 6):
            y_temp = old_fields.copy()
            for l in range(j):
                y_temp += dt * b[j][l] * k[l]
            k[j] = equations_FD(eq, y_temp, t + a[j]*dt, x_grid, y_grid, parameters, operators)

        # Combine stages to update field (5th order accurate)
        new_fields = old_fields.copy()
        for j in range(6):
            new_fields += dt * c[j] * k[j]

        fields = new_fields

        # Save history every t_rate steps
        if i % t_rate == 0:
            fields_history.append(fields.copy())
            time_grid.append(t)

    return fields, fields_history, time_grid


def RK4_ode(eq, vector, t_grid, parameters, dt, Nt):
    vector_history = []
    time_grid = []
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
            time_grid.append(t_grid[i])
    vector_history = np.array(vector_history)
    return vector, vector_history, time_grid


def RK4_FD_lyapunov(eq, fields, parameters, grids, dt, Nt, operators, t_rate, N_condit, frac_time):
    t_grid = grids[0]
    x_grid = grids[1]
    y_grid = grids[2]
    fields_history = []
    time_grid_lyap = []
    time_grid_fields = []
    lyap = []
    Nx = len(x_grid)
    I = np.eye(2 * Nx)
    U_init = np.random.rand(2 * Nx, N_condit) - 0.5
    Q, R = np.linalg.qr(U_init)
    dt_over_2 = 0.5 * dt

    for i in range(Nt - 1):
        old_fields = fields
        k_1 = equations_FD(eq, old_fields, t_grid[i], x_grid, y_grid, parameters, operators)
        k_2 = equations_FD(eq, old_fields + dt_over_2 * k_1, t_grid[i], x_grid, y_grid, parameters, operators)
        k_3 = equations_FD(eq, old_fields + dt_over_2 * k_2, t_grid[i], x_grid, y_grid, parameters, operators)
        k_4 = equations_FD(eq, old_fields + dt * k_3, t_grid[i], x_grid, y_grid, parameters, operators)
        new_fields = old_fields + dt * (k_1 + 2 * k_2 + 2 * k_3 + k_4) / 6
        fields = new_fields
        if i > (Nt - 1) - frac_time:
            J = jacobians_FD(eq, [new_fields[0], new_fields[1]], t_grid, x_grid, [0], parameters, operators)
            Q_new = time_propagator("I_Jdt_RK4", I, J, Q, dt)
            Q, R = np.linalg.qr(Q_new)
            Q = Q
            lyap.append(np.log(np.absolute(R.diagonal())) / dt)
            time_grid_lyap.append(t_grid[i])
        if i % t_rate == 0:
            fields_history.append(fields)
            time_grid_fields.append(t_grid[i])
    return time_grid_lyap, lyap, fields_history, time_grid_fields