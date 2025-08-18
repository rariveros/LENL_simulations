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

if __name__ == '__main__':
    eq = "PT_waveguide"
    t_rate = 1

    # PARAMETROS (GAMMA = POTENCIAL, K = ACOPLE, ALPHA = GAIN - LOSS)
    gamma = 10.0
    k = 1.0
    beta = 0.5
    delta = 0.0

    # PARAMETROS INICIALES (P = CANTIDAD CONSERVADA, X = PORCENTAJE DE INFORMACIÓN INICIAL EN DIMERO 2) ####### Vale pico esto, encuentra los puntos estacionarios como la gente
    Ns = np.array([16, 18, 20, 22])
    colors = ["r", "r", "g", "g"]
    theta_0 = 0.0
    n = 0
    phi_01 = np.arcsin(beta / k)
    phi_02 = np.pi - np.arcsin(beta / k)
    z_0 = [0.1]#np.arange(0, 1, 0.2)
    z = z_0
    phi = np.array([phi_01] * len(Ns))  # n * np.pi
    Nphi = len(phi)

    # CONDICIONES INICIALES EN TERMINOS DE P Y X
    U_init = []
    V_init = []
    COLOR = []
    for i in range(len(z)):
        for j in range(len(phi)):
            N = Ns[j]
            R1 = np.sqrt((N - z) / 2)
            R2 = np.sqrt((N + z) / 2)
            U_init.append(R1[i] * np.exp(1j * theta_0))
            V_init.append(R2[i] * np.exp(1j * (theta_0 + phi[j])))
            COLOR.append(colors[j])
    U_init = np.array(U_init)
    V_init = np.array(V_init)

    # DEFINICIÓN DE GRILLA TEMPORAL, ESPACIAL SE DEFINE COMO ARRAY CERO POR COMO FUNCIONA EL CODIGO
    [tmin, tmax, dt] = [0, 5000, 0.025]
    t_grid = np.arange(tmin, tmax + dt, dt)         # TEMPORAL GRID DEFINITION
    x_grid = np.array([0])           # SPATIAL GRID DEFINITION
    T = tmax
    Nt = t_grid.shape[0]

    # CONDICIONES INICIALES EN TERMINOS DE P Y X

    print(np.abs(U_init))
    print(np.abs(V_init))

    # EMPAQUETAMIENTO DE PARAMETROS PARA SIMULACIÓN
    operators = [0]
    fields_init = [U_init, V_init]
    grids = [t_grid, x_grid, 0]
    parameters_np = np.array([k, gamma, beta, delta])

    # SIMULACIÓN NUMERICA
    final_fields, fields_history, time_grid = RK4_FD(eq, fields_init, parameters_np, grids, dt, Nt, operators, t_rate)      #NUMERICAL SIMULATION

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
    arg_variable = np.angle(U_light * np.conjugate(V_light))

    ############### GENERIC INITIAL AND FINAL TIME ###############

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for i in range(len(N_variable[0, :])):
        ax.plot3D(N_variable[:, i], z_variable[:, i], np.unwrap(arg_variable, axis=0)[:, i], label='3D curve', color=COLOR[i], lw=0.2)
    ax.set_zlim([-np.pi, np.pi])
    ax.set_xlabel('$N_{1}$')
    ax.set_ylabel('$N_{2}$')
    ax.set_zlabel(r'$\phi$')
    plt.show()
    plt.close()