from back_process import *
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate
import os
import time
import datetime
from scipy.sparse import diags, hstack, vstack, identity, csr_matrix
from scipy.linalg import expm


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

    dir = "D:/mnustes_science/simulation_data/FD/network_chimeras/FIG01/erdos_renyi/k=0.0150/mean_degree=22.00"
    U_light = np.loadtxt(dir + '/U.txt', delimiter=',')
    V_light = np.loadtxt(dir + '/V.txt', delimiter=',')
    t_light = np.loadtxt(dir + '/T.txt', delimiter=',')
    x_grid = np.loadtxt(dir + '/X.txt', delimiter=',')
    laplacian = np.loadtxt(dir + '/L.txt', delimiter=',')
    parameters_np = np.loadtxt(dir + '/params.txt', delimiter=',')

    phase = np.arctan2(V_light, U_light)
    phase_light_wraped = phase
    phase_light = np.unwrap(phase_light_wraped, axis=0)
    module = np.sqrt(U_light ** 2 + V_light ** 2)

    dt = t_light[1] - t_light[0]
    Nt = len(t_light)
    Nx = len(x_grid)

    alpha = 0.4  # NONLINEAR COEFFICIENT
    mu = 0.1  # DISSIPATION
    gamma = 2.90  # 2.7  2.90                           # DRIVE STRENGTH
    k = 0.015  # 0.4216 #0.028                          # COUPLING
    w = 0.7  # FORCING FREQUENCY
    eq = 'duffing'  # OSCILLATOR TYPE
    t_rate = 1
    mean_degree = 22
    N_nodes = Nx

    [alpha, mu, gamma, k, w] = parameters_np
    L_csr = csr_matrix(laplacian)
    operators = [L_csr]

    n = 501
    p = mean_degree / (n - 1)

    t_init = 0
    t_final = 10000
    i_0 = np.argmin(np.abs(t_light - t_init))
    i_f = np.argmin(np.abs(t_light - t_final))
    power_threshold = 1.5

    file = dir
    k_str = f"{k:.{4}f}"

    if not os.path.exists(file):
        os.makedirs(file)

    #### LYAPUNOV SPECTRUM ####
    print("### LYAPUNOV ANALYISIS ###")
    now = datetime.datetime.now()
    print('Hora de Inicio: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
    time_init = time.time()

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

    x = np.arange(N_condit) / N_condit
    plt.plot(x, lyap_mean, c="k", linewidth=2)
    plt.scatter(x, lyap_mean, color="k")
    plt.hlines(0, 1, 0, color="k", linestyle="--")
    plt.title("Lyapunov Spectrum", size='15')
    plt.xlabel('Initial Conditions', size='15')
    plt.xticks(fontsize=15)
    plt.xlim(0, 1)
    plt.ylabel('$\\lambda$', size='20')
    plt.yticks(fontsize=15)
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig(file + '/lyap_espectrums.png', dpi=300)
    plt.close()

    arg_lyap_pos = lyap_mean >= 0.0
    Q_pos = (np.sum(np.abs(Q[:N_nodes, arg_lyap_pos]) + np.abs(Q[N_nodes:, arg_lyap_pos]), axis=1)) / (
        len(lyap_mean[arg_lyap_pos]))
    args_lyap_modes = np.argsort(Q_pos)

    plt.scatter(x_grid, Q_pos[args_lyap_modes], c="k", label="Sorted")
    plt.scatter(x_grid, Q_pos, c="r", label="Non-sorted")
    plt.xlabel("index", fontsize=20)
    plt.ylabel("Chaos Contribution", fontsize=20)
    plt.legend(loc="upper right")
    plt.xlim(0, N_nodes)
    plt.ylim(0, 1.1 * np.amax(Q_pos))
    plt.savefig(file + '/chaos_contribution.png', dpi=300)
    plt.close()

    now = datetime.datetime.now()
    print('Hora de Término: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
    time_fin = time.time()
    print(str(time_fin - time_init) + ' seg')

    np.savetxt(file + '/arg_lyap.txt', arg_lyap_pos, delimiter=',')
    np.savetxt(file + '/x_lyap.txt', x, delimiter=',')
    np.savetxt(file + '/lyap_spectrum.txt', lyap_mean, delimiter=',')

    #### SPATIOTEMPORAL DIAGRAMS ####

    fig, ((ax01, ax02), (ax03, ax04)) = plt.subplots(2, 2, figsize=(10, 8))
    cax_01 = ax01.pcolormesh(x_grid, t_light[i_0:i_f:100], U_light[i_0:i_f:100, :], cmap="turbo", shading='auto')
    cbar = fig.colorbar(cax_01)
    cbar.ax.tick_params(labelsize=13)
    cbar.set_label('$u(t)$', rotation=0, size=20, labelpad=-50, y=1.1)
    ax01.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
    ax01.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=False, labeltop=False, labelbottom=False)
    ax01.set_ylabel("$t$", fontsize=20)

    cax_02 = ax02.pcolormesh(x_grid, t_light[i_0:i_f:100], module[i_0:i_f:100, :], cmap="turbo", shading='auto')
    cbar = fig.colorbar(cax_02)
    cbar.ax.tick_params(labelsize=13)
    cbar.set_label('$r(t)$', rotation=0, size=20, labelpad=-50, y=1.1)
    ax02.tick_params(axis="y", direction="in", labelsize=15, left=False, right=True, labelleft=False, labelright=False)
    ax02.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=False, labeltop=False, labelbottom=False)

    cax_03 = ax03.pcolormesh(x_grid, t_light[i_0:i_f:100], U_light[i_0:i_f:100, args_lyap_modes], cmap="turbo", shading='auto')
    cbar = fig.colorbar(cax_03)
    cbar.ax.tick_params(labelsize=13)
    cbar.set_label('$u(t)$', rotation=0, size=20, labelpad=-50, y=1.1)
    ax03.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
    ax03.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=True)
    ax03.set_xlabel("$i$", fontsize=20)
    ax03.set_ylabel("$t$", fontsize=20)

    cax_04 = ax04.pcolormesh(x_grid, t_light[i_0:i_f:100], module[i_0:i_f:100, args_lyap_modes], cmap="turbo", shading='auto')
    cbar = fig.colorbar(cax_04)
    cbar.ax.tick_params(labelsize=13)
    cbar.set_label('$r(t)$', rotation=0, size=20, labelpad=-35, y=1.1)
    ax04.tick_params(axis="y", direction="in", labelsize=15, left=False, right=True, labelleft=False, labelright=False)
    ax04.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=True)
    ax04.set_xlabel("$i$", fontsize=20)

    plt.tight_layout()
    plt.savefig(file + '/spatiotemporals.png', dpi=200)
    plt.close()