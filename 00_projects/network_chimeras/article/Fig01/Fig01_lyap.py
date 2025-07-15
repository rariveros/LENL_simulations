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
    x_grid = np.loadtxt(dir + '/X.txt', delimiter=',')
    file = dir
    N_nodes = 501

    arg_lyap_pos = np.loadtxt(file + '/arg_lyap.txt', delimiter=',')
    x = np.loadtxt(file + '/x_lyap.txt', delimiter=',')
    lyap_mean = np.loadtxt(file + '/lyap_spectrum.txt', delimiter=',')
    Q_pos = np.loadtxt(file + '/Q_pos.txt', delimiter=',')

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

