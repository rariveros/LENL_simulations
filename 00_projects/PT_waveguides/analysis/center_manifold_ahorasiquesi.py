import matplotlib.pyplot as plt
import numpy as np

# DEFINICION DE INTEGRADOR TEMPORAL
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

# DEFINICION DE SISTEMA DE ECUACIONES
def equations_FD(eq, field_slices, t_i, x_grid, y_grid, parameters, operators): #ecuaciones
    if eq == "PT_waveguide":
        U = field_slices[0]
        V = field_slices[1]

        k = parameters[0]
        gamma = parameters[1]
        beta = parameters[2]
        delta = parameters[3]

        F = 1j * k * V - 1j * gamma * (U / (1 + np.abs(U) ** 2)) + (beta + 1j * delta) * U
        G = 1j * k * U - 1j * gamma * (V / (1 + np.abs(V) ** 2)) - (beta + 1j * delta) * V

        fields = np.array([F, G])
    return fields

def v0(N, beta, k):
    """Autovector teórico asociado al autovalor cero"""
    return np.array([1, 0, (2 * beta) / (N * np.sqrt(k**2 - beta**2))])

def EOM(k, gamma, beta, N, z, phi):
    dN = 2 * beta * z
    dz = 2 * k * np.sqrt(N ** 2 - z ** 2) * np.sin(phi) - 2 * beta * N
    dphi = - ((2 * k * z) / (np.sqrt(N ** 2 - z ** 2))) * np.cos(phi) + ((4 * gamma * z) / ((N + 2) ** 2 - z ** 2))
    return dN ** 2 + dz ** 2 + dphi ** 2

def center_manifold(gamma, k, beta, Ns, N_fix, phi_fix):
    phi_N = -(2 * beta / (np.sqrt(k ** 2 - beta ** 2))) * np.log(Ns / N_fix) + phi_fix
    return phi_N

def inverse_center_manifold(gamma, k, beta, phis, N_fix, phi_fix):
    N_phi = N_fix * np.exp((k / (2 * beta)) * (np.sin(phis) - np.sin(phi_fix)))#N_fix * np.exp(-(2 * beta / (np.sqrt(k ** 2 - beta ** 2))) ** (-1) * (phis - phi_fix))
    return N_phi

def jacobian(k, gamma, beta, N, z, phi):
    J11 = 0
    J12 = 2 * beta
    J13 = 0

    J21 = (2 * k * N / (N ** 2 - z ** 2) ** 0.5) * np.sin(phi) - 2 * beta
    J22 = - (2 * k * z / (N ** 2 - z ** 2) ** 0.5) * np.sin(phi)
    J23 = 2 * k * (N ** 2 - z ** 2) ** 0.5 * np.cos(phi)

    J31 = (2 * k * N * z / (N ** 2 - z ** 2) ** 1.5) * np.cos(phi) + 8 * gamma * z * (N + 2) / (
                (N + 2) ** 2 - z ** 2) ** 2
    J32 = -(2 * k * N ** 2 / (N ** 2 - z ** 2) ** 1.5) * np.cos(phi) + 4 * gamma * ((N + 2) ** 2 + z ** 2) / (
                (N + 2) ** 2 - z ** 2) ** 2
    J33 = (2 * k * z / np.sqrt(N ** 2 - z ** 2)) * np.sin(phi)

    J = np.array([
        [J11, J12, J13],
        [J21, J22, J23],
        [J31, J32, J33]
    ])
    return J

if __name__ == '__main__':
    gamma = 10.0
    k = 1.0
    dN = 0.2
    beta = 0.5
    delta = 0.0

    eq = "PT_waveguide"
    t_rate = 1

    Ns = np.arange(0.1, 50, dN)
    phis = np.arange(np.arcsin(beta / k), - np.pi - np.arcsin(beta / k), -0.01)
    N_fixed = np.array([14, 30, 45])
    phi_fixed = np.array([np.arcsin(beta / k)] * 3)

    f0 = center_manifold(gamma, k, beta, Ns, N_fixed[0], phi_fixed[0])
    f1 = center_manifold(gamma, k, beta, Ns, N_fixed[1], phi_fixed[1])
    f2 = center_manifold(gamma, k, beta, Ns, N_fixed[2], phi_fixed[2])

    g0 = inverse_center_manifold(gamma, k, beta, phis, N_fixed[0], phi_fixed[0])
    g1 = inverse_center_manifold(gamma, k, beta, phis, N_fixed[1], phi_fixed[1])
    g2 = inverse_center_manifold(gamma, k, beta, phis, N_fixed[2], phi_fixed[2])

    phi0 = center_manifold(gamma, k, beta, g0, N_fixed[0], phi_fixed[0])
    phi1 = center_manifold(gamma, k, beta, g1, N_fixed[1], phi_fixed[1])
    phi2 = center_manifold(gamma, k, beta, g2, N_fixed[2], phi_fixed[2])
    phases = [- np.pi - np.arcsin(beta / k), np.arcsin(beta / k), - np.arcsin(beta / k) + np.pi, np.arcsin(beta / k) + 2 * np.pi]
    colors_a = []
    colors_b = []
    #fig, ax = plt.subplots(5, 1, figsize=(4, 8), dpi=100)
    for i in range(len(Ns)):
        jacobian_ai = jacobian(k, gamma, beta, Ns[i], 0, phases[0])
        jacobian_bi = jacobian(k, gamma, beta, Ns[i], 0, phases[1])
        jacobian_ci = jacobian(k, gamma, beta, Ns[i], 0, phases[2])
        jacobian_di = jacobian(k, gamma, beta, Ns[i], 0, phases[3])
        eigvalues_a, eigvectors_a = np.linalg.eig(jacobian_ai)
        eigvalues_b, eigvectors_b = np.linalg.eig(jacobian_bi)
        eigvalues_c, eigvectors_c = np.linalg.eig(jacobian_ci)
        eigvalues_d, eigvectors_d = np.linalg.eig(jacobian_di)
        EOM_01 = EOM(k, gamma, beta, Ns[i], 0,  phases[0])
        EOM_02 = EOM(k, gamma, beta, Ns[i], 0,  phases[1])
        EOM_03 = EOM(k, gamma, beta, Ns[i], 0,  phases[2])
        EOM_04 = EOM(k, gamma, beta, Ns[i], 0,  phases[3])

        if np.sum(np.abs(np.real(eigvalues_a))) > 0.01:
            colors_a.append(1)
        else:
            colors_a.append(0)
        if np.sum(np.abs(np.real(eigvalues_b))) > 0.01:
            colors_b.append(1)
        else:
            colors_b.append(0)
        #ax[1].scatter([Ns[i]] * len(eigvalues_a), np.real(eigvalues_a), c="b")
        #ax[1].scatter([Ns[i]] * len(eigvalues_a), np.imag(eigvalues_a), c="r")
        #ax[2].scatter([Ns[i]] * len(eigvalues_b), np.real(eigvalues_b), c="b")
        #ax[2].scatter([Ns[i]] * len(eigvalues_b), np.imag(eigvalues_b), c="r")
        #ax[3].scatter([Ns[i]] * len(eigvalues_c), np.real(eigvalues_c), c="b")
        #ax[3].scatter([Ns[i]] * len(eigvalues_c), np.imag(eigvalues_c), c="r")
        #ax[4].scatter([Ns[i]] * len(eigvalues_d), np.real(eigvalues_d), c="b")
        #ax[4].scatter([Ns[i]] * len(eigvalues_d), np.imag(eigvalues_d), c="r")
        #ax[0].scatter(Ns[i], EOM_01, c="k")
        #ax[0].scatter(Ns[i], EOM_02, c="k")
        #ax[0].scatter(Ns[i], EOM_03, c="k")
        #ax[0].scatter(Ns[i], EOM_04, c="k")
    #plt.show()
    #plt.close()
    ls_a = np.argwhere(np.diff(np.abs(np.array(colors_a))) > 0.001)
    ls_b = np.argwhere(np.abs(np.diff(np.array(colors_b))) > 0.001)
    print(ls_a)
    print(ls_b)
    each = 6
    g0 = g0[::each]
    g1 = g1[::each]
    g2 = g2[::each]
    phi0 = phi0[::each]
    phi1 = phi1[::each]
    phi2 = phi2[::each]
    phis = phis[::each]


    N_max = 50
    fig, ax = plt.subplots(1, 1, figsize=(5, 4), dpi=150)
    ax.plot(Ns, f0, color="r", zorder=5)
    ax.plot(Ns, f1, color="g", zorder=5)
    ax.plot(Ns, f2, color="b", zorder=5)

    ax.hlines(phases[0], 0, N_max, colors="k", linestyles="-")
    ax.hlines(phases[1], 0, Ns[98], colors="k", linestyles="--")
    ax.hlines(phases[1], Ns[98], N_max, colors="k", linestyles="-")
    ax.hlines(phases[2], 0, N_max, colors="k", linestyles="-")
    ax.hlines(phases[3], 0, Ns[98], colors="k", linestyles="--")
    ax.hlines(phases[3], Ns[98], N_max, colors="k", linestyles="-")

    ax.set_xlim(0, N_max)
    ax.set_ylim(0, np.pi)
    ax.set_xlabel("$N_0$")
    ax.set_ylabel("$\phi$")
    plt.show()
    plt.close()

    #plt.plot(g2, phi2, color="k")
    #plt.scatter(g2[0], phi2[0], c="r")
    Ns = Ns[80::each]
    U_init = []
    V_init = []
    for i in range(len(Ns)):
        N_init = Ns[i] + 1
        z_init = 0
        phi_init = np.arcsin(beta/k)#phis[i]
        R1_init = np.sqrt((N_init - z_init) / 2)
        R2_init = np.sqrt((N_init + z_init) / 2)
        U_init.append(R1_init)
        V_init.append(R2_init * np.exp(1j * (phi_init)))
    U_init = np.array(U_init)
    V_init = np.array(V_init)

    # DEFINICIÓN DE GRILLA TEMPORAL, ESPACIAL SE DEFINE COMO ARRAY CERO POR COMO FUNCIONA EL CODIGO
    [tmin, tmax, dt] = [0, 200, 0.025]
    t_grid = np.arange(tmin, tmax + dt, dt)  # TEMPORAL GRID DEFINITION
    x_grid = np.array([0])  # SPATIAL GRID DEFINITION
    T = tmax
    Nt = t_grid.shape[0]

    # EMPAQUETAMIENTO DE PARAMETROS PARA SIMULACIÓN
    operators = [0]
    fields_init = [U_init, V_init]
    grids = [t_grid, x_grid, 0]
    parameters_np = np.array([k, gamma, beta, delta])

    # SIMULACIÓN NUMERICA
    final_fields, fields_history, time_grid = RK4_FD(eq, fields_init, parameters_np, grids, dt, Nt, operators, t_rate)

    # REOBTENIENDO DATOS DE SIMULACIÓN
    U = np.array(fields_history)[:, 0]
    V = np.array(fields_history)[:, 1]
    lightness = 1
    U_light = U[0::lightness]
    V_light = V[0::lightness]
    t_light = np.array(time_grid[0::lightness])
    P1 = np.abs(U_light) ** 2
    P2 = np.abs(V_light) ** 2

    N_variable = P1 + P2
    z_variable = P2 - P1
    arg_variable = -np.angle(U_light * np.conjugate(V_light))

    fig = plt.figure(dpi=160)
    ax = fig.add_subplot(111, projection='3d')
    for i in range(len(N_variable[0, :])):
        ax.plot3D(N_variable[:, i], z_variable[:, i], -np.unwrap(arg_variable[:, i]), label='3D curve', color="b", lw=2, alpha=0.8)
        # ax.scatter3D(N_array[i], 0, phi_array[i], c="r")

    ax.set_zlim([-0 * np.pi, -2 * np.pi])
    ax.set_xlabel('$N$')
    ax.set_ylabel('$z$')
    ax.set_zlabel(r'$-\phi$')
    ax.view_init(elev=-20, azim=123, roll=-176)
    plt.show()
    # plt.savefig("alo.png", dpi=300)
    plt.close()