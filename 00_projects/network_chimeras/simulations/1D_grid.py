import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sparse
from scipy.sparse import diags, hstack, vstack, identity
from scipy.linalg import expm

def RK4_FD(eq, fields, parameters, grids, dt, Nt, operators, t_rate): #implementa rouge-kutta
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

def sparse_DD_neumann(Nx, dx):
    data = np.ones((3, Nx))
    data[1] = -2 * data[1]
    diags = [-1, 0, 1]
    D2 = sparse.spdiags(data, diags, Nx, Nx) / (dx ** 2)
    D2 = sparse.lil_matrix(D2)
    # Condiciones de borde de Neumann: ajustar la primera y última fila
    D2[0, 0] = -1 / (dx ** 2)
    D2[0, 1] = 1 / (dx ** 2)
    D2[-1, -1] = -1 / (dx ** 2)
    D2[-1, -2] = 1 / (dx ** 2)
    return D2.tocsr()

def sparse_DD(Nx, dx):
    data = np.ones((3, Nx))
    data[1] = -2 * data[1]
    diags = [-1, 0, 1]
    D2 = sparse.spdiags(data, diags, Nx, Nx) / (dx ** 2)
    D2 = sparse.lil_matrix(D2)
    D2[0, -1] = 1 / (dx ** 2)
    D2[-1, 0] = 1 / (dx ** 2)
    return D2

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
    alpha = 0.4                                     # NONLINEAR COEFFICIENT
    mu = 0.1                                        # DISSIPATION
    gamma = 2.75#2.7  2.90                                  # DRIVE STRENGTH                                # COUPLING (0.42)
    w = 0.7
    eq = 'duffing'
    t_rate = 1

    # Grid definition
    N_nodes = 300
    [tmin, tmax, dt] = [0, 3000, 0.05]
    t_grid = np.arange(tmin, tmax + dt, dt)  # TEMPORAL GRID DEFINITION
    [xmin, xmax, dx] = [0, N_nodes, 1]
    x_grid = np.arange(xmin, xmax, dx)  # SPATIAL GRID DEFINITION

    sigma = 30
    k = 0.3 * (2 * np.exp(-(x_grid - 150) ** 2 / (2 * sigma ** 2)) + 1)  # 0.4216 #0.028

    T = tmax
    Nt = t_grid.shape[0]
    Nx = x_grid.shape[0]

    # Initial Conditions
    U_init = 0.1 * np.ones(Nx) + 0.025 * np.random.rand(Nx)
    initial_quimera = 100
    arg_chimera = [150] #np.arange(146, 155)  # INITIAL QUIMERA INDEX
    for i in arg_chimera:
        U_init[i] = 0.1
    V_init = 0.025 * np.random.rand(Nx)

    # Empaquetamiento de parametros, campos y derivadas para integración
    laplacian_matrix = sparse_DD(Nx, dx)
    operators = [laplacian_matrix]
    fields_init = [U_init, V_init]
    grids = [t_grid, x_grid, 0]
    parameters_np = [alpha, mu, gamma, k, w]

    final_fields, fields_history, time_grid = RK4_FD(eq, fields_init, parameters_np, grids, dt, Nt, operators, t_rate)

    # Reobteniendo campos y variables
    U = np.array(fields_history)[:, 0]
    V = np.array(fields_history)[:, 1]
    phase = np.arctan2(V, U)  # Harmonic oscillator can be writen as R * exp (i * phase) with R = SQRT(U ** 2 + V ** 2) and phase = arctan(V/U)
    lightness = 1
    U_light = U[0::lightness]
    V_light = V[0::lightness]
    phase_light_wraped = phase[0::lightness]
    phase_light = np.unwrap(phase_light_wraped, axis=0)[0::lightness]
    t_light = np.array(time_grid[0::lightness])
    module = np.sqrt(U_light ** 2 + V_light ** 2)  ## AMPLITUDE / MODULE

    print("simulation done")
    plt.plot(x_grid, k, color="k")
    plt.show()
    plt.close()

    t_init = 2800
    t_final = 3000
    i_0 = np.argmin(np.abs(t_light - t_init))
    i_f = np.argmin(np.abs(t_light - t_final))

    fig, (ax01, ax02) = plt.subplots(2, 1, figsize=(14, 12))
    cax_01 = ax01.pcolormesh(x_grid, t_light[i_0:i_f], U_light[i_0:i_f, :], cmap="turbo", shading='auto')
    cbar = fig.colorbar(cax_01)
    cbar.ax.tick_params(labelsize=13)
    cbar.set_label('$u(t)$', rotation=0, size=20, labelpad=-50, y=1.1)
    ax01.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
    ax01.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=True)
    ax01.set_xlabel("$i$", fontsize=20)
    ax01.set_ylabel("$t$", fontsize=20)

    cax_02 = ax02.pcolormesh(x_grid, t_light[i_0:i_f], module[i_0:i_f, :], cmap="turbo", shading='auto')
    cbar = fig.colorbar(cax_02)
    cbar.ax.tick_params(labelsize=13)
    cbar.set_label('$r(t)$', rotation=0, size=20, labelpad=-50, y=1.1)
    ax02.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
    ax02.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=True)
    ax02.set_xlabel("$i$", fontsize=20)
    ax02.set_ylabel("$t$", fontsize=20)

    figure_dir = "C:/Users/Rafa/PycharmProjects/LENL_simulations/00_projects/network_chimeras/simulations/figures"
    plt.savefig(figure_dir + '/fig_chimeras.png', dpi=300)

    from scipy import integrate

    ############## AVERAGES AND CHARACTERIZATION OF THE STATIONARY FINAL STATE N = (10000, 20000)##############

    average_phase = integrate.simpson(phase_light_wraped[i_0:i_f], t_grid[i_0:i_f], axis=0) / (
                t_grid[i_f] - t_grid[i_0])  ## AVERAGE WRAPPED PHASE
    omega_average = integrate.simpson(np.diff(phase_light, axis=0)[i_0:i_f] / dt, t_grid[i_0:i_f], axis=0) / (
                t_grid[i_f] - t_grid[i_0])  ## AVERAGE FREQUENCY
    sigma2_phase = integrate.simpson((phase_light[i_0:i_f] - average_phase) ** 2, t_grid[i_0:i_f], axis=0) / (
                t_grid[i_f] - t_grid[i_0])  ## FREQUENCY STANDARD DEVIATION
    sigma2 = integrate.simpson(((np.diff(phase_light, axis=0)[i_0:i_f] / dt) - omega_average) ** 2, t_grid[i_0:i_f],
                               axis=0) / (t_grid[i_f] - t_grid[i_0])  ## FREQUENCY STANDARD DEVIATION
    average_module = integrate.simpson(module[i_0:i_f], t_grid[i_0:i_f], axis=0) / (t_grid[i_f] - t_grid[i_0])

    args = np.argsort(average_module)  ## ACA SE ORDENA

    fig, (ax01, ax02, ax03, ax04, ax05) = plt.subplots(5, 1, figsize=(10, 10))

    ax01.scatter(x_grid, average_phase, color="k")
    ax01.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
    ax01.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=False)
    ax01.set_ylabel("$\\langle \\theta_i \\rangle$", fontsize=20)

    ax02.scatter(x_grid, average_module, color="k")
    ax02.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
    ax02.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=False)
    ax02.set_ylabel("$\\langle  r_i \\rangle$", fontsize=20)

    ax03.scatter(x_grid, omega_average, color="k")
    ax03.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
    ax03.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=False)
    ax03.set_ylabel("$\\langle \\omega_i \\rangle$", fontsize=20)

    ax04.scatter(x_grid, np.sqrt(sigma2_phase), color="k")
    ax04.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
    ax04.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=False)
    ax04.set_ylabel("$\\sigma_i$", fontsize=20)

    ax05.scatter(x_grid, phase_light_wraped[-1, :], color="k")
    ax05.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
    ax05.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=True)
    ax05.set_ylabel("$\\theta_i(t_f)$", fontsize=20)
    ax05.set_xlabel("$i$", fontsize=20)

    plt.show()

    import scipy.linalg

    T_00 = int(0.98 * Nt)
    T_01 = int(1.00 * Nt)
    I = np.eye(2 * Nx)

    N_condit = 2 * Nx
    U_init = np.random.rand(N_condit, N_condit) - 0.5
    Q, _ = np.linalg.qr(U_init)  # Only Q is needed

    # Preallocate Lyapunov exponent storage for speed
    lyap_01 = np.zeros((T_01 - T_00 - 1, N_condit))

    for i in range(T_00, T_01 - 1):
        U_1, U_2 = U_light[i, :], V_light[i, :]
        J = jacobians_FD(eq, [U_1, U_2], t_grid, x_grid, [0], parameters_np, operators)

        # Compute the time propagator and update Q
        Q_new = time_propagator("I_Jdt_RK4", I, J, Q, dt)
        Q, R = scipy.linalg.qr(Q_new, overwrite_a=True, mode='economic')  # QR decomposition

        # Store log of absolute diagonal of R divided by dt
        lyap_01[i - T_00, :] = np.log(np.abs(np.diag(R))) / dt

    # Compute the mean Lyapunov exponents
    lyap_mean_01 = np.mean(lyap_01, axis=0)

    x = np.arange(N_condit) / N_condit
    plt.plot(x, lyap_mean_01, c="k", linewidth=2)
    plt.scatter(x, lyap_mean_01, color="k")
    plt.hlines(0, 1, 0, color="k", linestyle="--")
    plt.title("Lyapunov Spectrum", size='15')
    plt.xlabel('Initial Conditions', size='15')
    plt.xticks(fontsize=15)
    plt.ylabel('$\\lambda$', size='20')
    plt.yticks(fontsize=15)
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig('lyap_espectrums.png', dpi=300)
    plt.show()
    plt.close()

    arg_lyap_pos = lyap_mean_01 >= 0.0
    Q_pos = (np.sum(np.abs(Q[:N_nodes, arg_lyap_pos]) + np.abs(Q[N_nodes:, arg_lyap_pos]), axis=1)) / (
        len(lyap_mean_01[arg_lyap_pos]))
    args_lyap_modes = np.argsort(Q_pos)

    plt.plot(t_grid[T_00:T_01 - 1], lyap_01[:, 30])
    plt.show()
    plt.close()

    plt.scatter(x_grid, Q_pos[args_lyap_modes], c="k", label="Sorted")
    plt.scatter(x_grid, Q_pos, c="r", label="Non-sorted")
    plt.xlabel("index", fontsize=20)
    plt.ylabel("Chaos Contribution", fontsize=20)
    plt.legend(loc="upper right")
    plt.xlim(0, N_nodes)
    plt.ylim(0, 1.1 * np.amax(Q_pos))
    plt.show()
    plt.close()