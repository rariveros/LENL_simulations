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

if __name__ == '__main__':
    gamma = 10.0
    k = 1.0
    dN = -0.001
    Ns = np.arange(20, 0.01, dN)
    beta = 0.5
    delta = 0.0

    eq = "PT_waveguide"
    t_rate = 1

    # Inicializa listas para guardar los datos
    N_list = []
    z_list = []
    phi_list = []
    COLOR = []

    for i in range(len(Ns)):
        if i == 0:
            N = Ns[i]
            phi = np.arcsin(beta / k)
        else:
            N = N + dN
            phi = phi + dN * (2 * beta) / (N * np.sqrt(k ** 2 - beta ** 2))

        a = 2 * beta
        c = 2 * N * np.sqrt(k ** 2 - beta ** 2)
        d = 4 * gamma / (N + 2) ** 2 - 2 * np.sqrt(k ** 2 - beta ** 2) / N

        color_01 = "r" if np.real(np.sqrt(c * d + 0j)) > 0.01 else "b"
        color_02 = "r" if np.real(np.sqrt(c * d - 2 * a ** 2 + 0j)) > 0.01 else "b"
        COLOR.append(color_02)

        z = 0  # z siempre cero
        N_list.append(N)
        z_list.append(z)
        phi_list.append(phi)

    # Convertir a arrays
    N_array = np.array(N_list)
    z_array = np.array(z_list)
    phi_array = np.array(phi_list)

    # Derivadas por diferencias finitas
    dN_vec = np.diff(N_array)
    dz_vec = np.diff(z_array)
    dphi_vec = np.diff(phi_array)

    dxds_num = np.vstack([dN_vec, dz_vec, dphi_vec]).T / dN

    # Autovector v0 evaluado en cada punto (sin el último para alinear con derivadas)
    v0_array = np.array([v0(N, beta, k) for N in N_array[:-1]])

    # Error entre derivada y autovector
    errors = np.linalg.norm(dxds_num - v0_array, axis=1)

    # Plot del error
    plt.plot(N_array[:-1], errors)
    plt.xlabel("$N$")
    plt.ylabel("Error $||d\mathbf{x}/ds - \mathbf{v}_0||$")
    plt.title("Error entre vector tangente y autovector nulo")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    #ax.set_xlabel('$N$')
    #ax.set_ylabel('$z$')
    #ax.set_zlabel(r'$\phi$')
    #plt.show()
    #plt.close()

    each = 500
    # Convierte a arrays de numpy
    N_array = np.array(N_list[::each])
    z_array = np.array(z_list[::each])
    phi_array = np.array(phi_list[::each])
    COLOR = COLOR[::each]

    U_init = []
    V_init = []
    for i in range(len(N_array)):
        N_init = N_array[i]
        z_init = z_array[i]
        phi_init = phi_array[i]
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

    # CONDICIONES INICIALES EN TERMINOS DE P Y X

    print(np.abs(U_init))
    print(np.abs(V_init))

    # EMPAQUETAMIENTO DE PARAMETROS PARA SIMULACIÓN
    operators = [0]
    fields_init = [U_init, V_init]
    grids = [t_grid, x_grid, 0]
    parameters_np = np.array([k, gamma, beta, delta])

    # SIMULACIÓN NUMERICA
    final_fields, fields_history, time_grid = RK4_FD(eq, fields_init, parameters_np, grids, dt, Nt, operators, t_rate)  # NUMERICAL SIMULATION

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

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for i in range(len(N_variable[0, :])):
        ax.plot3D(N_variable[:, i], z_variable[:, i], np.unwrap(arg_variable[:, i]), label='3D curve', color="k", lw=0.2)
        #ax.scatter3D(N_array[i], 0, phi_array[i], c="r")

    ax.set_zlim([-np.pi, np.pi])
    ax.set_xlabel('$N$')
    ax.set_ylabel('$z$')
    ax.set_zlabel(r'$\phi$')
    plt.show()
    #plt.savefig("alo.png", dpi=300)
    plt.close()
