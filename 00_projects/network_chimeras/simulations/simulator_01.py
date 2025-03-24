import os
import time
import datetime
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
import scipy.sparse as sparse

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


def equations_FD(eq, field_slices, t_i, x_grid, y_grid, parameters, operators):
    if eq == 'duffing':
        U = field_slices[0]
        V = field_slices[1]

        alpha = parameters[0]
        mu = parameters[1]
        gamma = parameters[2]
        k = parameters[3]
        w = parameters[4]
        DD = operators[0]

        ddU = Der(DD, U)

        F = V
        G = - U + alpha * U ** 3 - U ** 5 - mu * V + gamma * np.cos(w * t_i) + k * ddU

        fields = np.array([F, G])
    return fields


def Der(D, f):
    d_f = D @ f
    return d_f

def sparse_DD(Nx, dx):
    data = np.ones((3, Nx))
    data[1] = -2 * data[1]
    diags = [-1, 0, 1]
    D2 = sparse.spdiags(data, diags, Nx, Nx) / (dx ** 2)
    D2 = sparse.lil_matrix(D2)
    D2[0, -1] = 1 / (dx ** 2)
    D2[-1, 0] = 1 / (dx ** 2)
    return D2

def sparse_D_neumann(Nx, dx):
    data = np.ones((2, Nx))
    data[0] = -1 * data[0]
    diags = [0, 1]
    D1 = sparse.spdiags(data, diags, Nx, Nx) / dx
    D1 = sparse.lil_matrix(D1)
    D1[0, 0] = 0
    D1[-1, -1] = 0
    return D1


if __name__ == '__main__':

    # Definiendo parámetros
    project_name = '/duffing'
    disc = 'D:/'
    route = 'mnustes_science/simulation_data/FD'    # SAVE DIRECTORY
    eq = 'duffing'                                  # EQUATION OF MOTION
    t_rate = 10                                     # SAVERATE IN SIMULATION

    alpha = 0.4                                     # NONLINEAR COEFFICIENT
    mu = 0.1                                        # DISSIPATION
    gamma = 2.75                                    # DRIVE STRENGTH
    k = 0.42                                        # COUPLING
    w = 0.7                                         # DRIVE FREQ

    ti = 0.

    # Grid definition
    [tmin, tmax, dt] = [0, 300, 0.01]
    t_grid = np.arange(tmin, tmax + dt, dt)         # TEMPORAL GRID DEFINITION
    [xmin, xmax, dx] = [0, 300, 1]
    x_grid = np.arange(xmin, xmax, dx)              # SPATIAL FRID DEFINITION

    T = tmax
    Nt = t_grid.shape[0]
    Nx = x_grid.shape[0]

    # Initial Conditions
    U_init = 0.1 * np.ones(Nx)
    d = 80                                          # INITIAL QUIMERA SIZE
    for i in range(Nx):
        if i > int(150 - d / 2) and i < int(150 + d / 2):
            U_init[i] = 2
    U_init = U_init + 0.01 * (np.random.rand(Nx) - 0.5)
    V_init = 0 * np.random.rand(Nx)

    DD = sparse_DD(Nx, dx)                          # DEFINITION OF COUPLING MATRIX
    operators = np.array([DD])

    # Empaquetamiento de parametros, campos y derivadas para integración
    fields_init = [U_init, V_init]
    grids = [t_grid, x_grid, 0]
    parameters_np = np.array([alpha, mu, gamma, k, w])

    # SIMULATION
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
    subfile = "/test"                               # SUBDIRECTORY NAME
    if not os.path.exists(file + subfile):
        os.makedirs(file + subfile)

    np.savetxt(file + subfile + '/U.txt', U, delimiter=',')
    np.savetxt(file + subfile + '/V.txt', V, delimiter=',')
    np.savetxt(file + subfile + '/parameters.txt', parameters_np, delimiter=',')
    np.savetxt(file + subfile + '/T.txt', t_grid, delimiter=',')

    phase = np.angle(signal.hilbert(U, axis=0))
    phase = np.unwrap(phase, axis=0)
    lightness = 1
    U_light = U[0::lightness]
    V_light = V[0::lightness]
    phase_light = phase[0::lightness]
    t_light = time_grid[0::lightness]
    R = (1 / Nx) * np.abs(np.sum(np.exp(1j * phase), axis=1))

    ########## PLOTS ##########
    t_init = 250
    t_final = 300


    fig, (ax01, ax02) = plt.subplots(2, 1, figsize=(4, 6))
    cax_01 = ax01.pcolormesh(x_grid, t_light, U_light, cmap="turbo", shading='auto')
    cbar = fig.colorbar(cax_01)
    cbar.ax.tick_params(labelsize=13)
    ax01.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
    ax01.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=False)
    ax01.set_ylabel("$t$", fontsize=20)
    ax01.set_ylim(t_init, t_final)

    cax_02 = ax02.pcolormesh(x_grid, t_light[:-1], np.diff(phase_light, axis=0) / dt, cmap="turbo", shading='auto')
    cbar = fig.colorbar(cax_02)
    cbar.ax.tick_params(labelsize=13)
    ax02.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
    ax02.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=True)
    ax02.set_ylabel("$t$", fontsize=20)
    ax02.set_xlabel("$i$", fontsize=20)
    ax02.set_ylim(t_init, t_final)

    plt.savefig(file + subfile + "/test_01.png", dpi=300, bbox_inches='tight')
    plt.close()

    fig, (ax01, ax02, ax03) = plt.subplots(3, 1, figsize=(4, 4))
    ax01.plot(t_light, U_light[:, 0], color="k")
    ax01.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
    ax01.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=False)
    ax01.set_ylabel("$x_0$", fontsize=20)
    ax01.set_xlim(t_init, t_final)

    ax02.plot(t_light, U_light[:, 150], color="k")
    ax02.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
    ax02.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=False)
    ax02.set_ylabel("$x_{150}$", fontsize=20)
    ax02.set_xlim(t_init, t_final)

    ax03.plot(t_light, R, color="k")
    ax03.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
    ax03.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=True)
    ax03.set_xlabel("$t$", fontsize=20)
    ax03.set_ylabel("$R(t)$", fontsize=20)
    ax03.set_xlim(t_init, t_final)

    plt.savefig(file + subfile + "/test_02.png", dpi=300, bbox_inches='tight')
    plt.close()

    fig, (ax01, ax02, ax03) = plt.subplots(3, 1, figsize=(4, 4))
    ax01.plot(t_light, phase_light[:, 0], color="k")
    ax01.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
    ax01.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=False)
    ax01.set_ylabel("$x_0$", fontsize=20)
    ax01.set_xlim(0, 50)

    ax02.plot(t_light, phase_light[:, 150], color="k")
    ax02.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
    ax02.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=False)
    ax02.set_ylabel("$x_{150}$", fontsize=20)
    ax02.set_xlim(0, 50)

    ax03.plot(t_light, R, color="k")
    ax03.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
    ax03.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=True)
    ax03.set_xlabel("$t$", fontsize=20)
    ax03.set_ylabel("$R(t)$", fontsize=20)
    ax03.set_xlim(0, 50)

    plt.savefig(file + subfile + "/test_03.png", dpi=300, bbox_inches='tight')
    plt.close()