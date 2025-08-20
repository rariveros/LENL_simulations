from functions import *
from back_process import *
from matplotlib.lines import lineStyles

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

# Define la función f_n
def f_n(z, n, N, k, gamma, delta):
    term1 = (-1) ** n * k * z * ((N + 2)**2 - z**2)
    term2 = -2 * gamma * z * np.sqrt(N**2 - z**2)
    term3 = delta * ((N + 2)**2 - z**2) * np.sqrt(N**2 - z**2)
    return term1 + term2 + term3

def find_roots_vs_N(Ns, k, gamma, delta, dz):
    ROOTS_0 = []
    ROOTS_1 = []

    for N in Ns:
        z_vals = np.arange(-N + dz, N - dz, dz)

        for n, ROOTS in zip([0, 1], [ROOTS_0, ROOTS_1]):
            f_vals = [f_n(z, n, N, k, gamma, delta) for z in z_vals]
            f_vals = np.array(f_vals)

            roots = []
            sign_ref = np.sign(f_vals[0])

            for i in range(1, len(f_vals)):
                sign_now = np.sign(f_vals[i])
                if sign_now != sign_ref and sign_now != 0:
                    z_root = 0.5 * (z_vals[i-1] + z_vals[i])
                    roots.append(z_root)
                    sign_ref = sign_now

            ROOTS.append(roots)

    # Aplanar resultados
    N_vals_0, z_roots_0 = [], []
    N_vals_1, z_roots_1 = [], []

    for i, root_list in enumerate(ROOTS_0):
        for r in root_list:
            N_vals_0.append(Ns[i])
            z_roots_0.append(r)

    for i, root_list in enumerate(ROOTS_1):
        for r in root_list:
            N_vals_1.append(Ns[i])
            z_roots_1.append(r)

    return (np.array(N_vals_0), np.array(z_roots_0),
            np.array(N_vals_1), np.array(z_roots_1))

if __name__ == '__main__':
    eq = "PT_waveguide"
    t_rate = 1.0

    # PARAMETROS (GAMMA = POTENCIAL, K = ACOPLE, ALPHA = GAIN - LOSS)
    gamma = 10
    k = 1.0
    beta = 0.0
    delta = 0.0

    # PARAMETROS INICIALES (P = CANTIDAD CONSERVADA, X = PORCENTAJE DE INFORMACIÓN INICIAL EN DIMERO 2) ####### Vale pico esto, encuentra los puntos estacionarios como la gente
    Ns = [26.0] #[5.0, 10.0, 15.0, 20.0, 25.0, 30.0]
    theta_0 = 0.0

    [tmin, tmax, dt] = [0, 300, 0.025]
    t_grid = np.arange(tmin, tmax + dt, dt)         # TEMPORAL GRID DEFINITION
    x_grid = np.array([0])           # SPATIAL GRID DEFINITION
    T = tmax
    Nt = t_grid.shape[0]

    n = 0
    z0 = 0.0 #
    U_init = []
    V_init = []
    phi_initial = []
    phi = np.arange(-2 * np.pi, 2 * np.pi, 0.025 * np.pi) #[-2 * np.pi, np.pi, 0.0, -np.pi, 2 * np.pi] #np.arange(-2 * np.pi, 2.01 * np.pi, 0.5 * np.pi)
    for i in range(len(Ns)):
        N = Ns[i]
        z_plus = ((N + 2) ** 2 - (2 * gamma ** 2) / (k**2) + ((2 * gamma * (-1) ** n) / (k)) * (gamma ** 2 / (k ** 2) - 4 * (1 + N)) ** 0.5) ** 0.5
        z_minus = ((N + 2) ** 2 - (2 * gamma ** 2) / (k**2) - ((2 * gamma * (-1) ** n) / (k)) * (gamma ** 2 / (k ** 2) - 4 * (1 + N)) ** 0.5) ** 0.5
        dz = N / 12
        z = np.array([0])#np.arange(-N + dz, N - dz, dz) #
         #[-0.1, 0.1] #np.arange(0, np.pi, np.pi * 0.99) #n * np.pi
        Nphi = len(phi)
        R0 = np.sqrt(N / 2)
        R1 = np.sqrt((N - z) / 2)
        R2 = np.sqrt((N + z) / 2)
        for i in range(len(phi)):
            U_init.append(R1 * np.exp(1j * theta_0))
            V_init.append(R2 * np.exp(1j * (theta_0 + phi[i])))
            if np.abs(phi[i]) == 2 * np.pi:
                phi_initial.append(phi[i] * np.ones(len(z)))
            else:
                phi_initial.append(0 * np.ones(len(z)))
    # DEFINICIÓN DE GRILLA TEMPORAL, ESPACIAL SE DEFINE COMO ARRAY CERO POR COMO FUNCIONA EL CODIGO

    # CONDICIONES INICIALES EN TERMINOS DE P Y X
        #U_init.append(R1[::1] * np.exp(1j * theta_0))
        #V_init.append(R2[::1] * np.exp(1j * theta_0))
        #U_init.append(R1 * np.exp(1j * theta_0))
        #V_init.append(R2 * np.exp(1j * (theta_0 + np.pi)))
    #U_init.append(R0 * np.exp(1j * theta_0 + np.zeros(Nphi)))
    #V_init.append(R0 * np.exp(1j * (theta_0 + phi)))


    U_init = np.array([x for arr in U_init for x in arr])
    V_init = np.array([x for arr in V_init for x in arr])
    phi_initial = np.array([x for arr in phi_initial for x in arr])

    #z_extra = [2, 4, 6, 8, 10, 11, 12, 13, 14]
    #phi_extra = [0, 0, 0, 0, 0, 0, 0, 0, 0]

    z_extra = []
    phi_extra = []
    for i in range(len(z_extra)):
        R1 = np.sqrt((N - z_extra[i]) / 2)
        R2 = np.sqrt((N + z_extra[i]) / 2)
        U_init = np.append(U_init, R1 * np.exp(1j * theta_0))
        V_init = np.append(V_init, R2 * np.exp(1j * (theta_0 + phi_extra[i])))
        phi_initial = np.append(phi_initial, phi_extra[i])

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
    arg_variable = np.unwrap(arg_variable, axis=0)
    for i in range(len(arg_variable[0, :])):
        arg_variable[:, i] = arg_variable[:, i] + phi_initial[i]

    lower = -np.pi
    upper = np.pi
    mask = (arg_variable < lower) | (arg_variable > upper)

    #N_variable[mask] = np.nan
    #z_variable[mask] = np.nan
    #arg_variable[mask] = np.nan

    n = 0
    z_pp_0 = +((N + 2) ** 2 - (2 * gamma ** 2) / (k ** 2) + ((2 * gamma) / (k ** 2)) * (
                gamma ** 2 - 4 * k * (-1) ** n * (1 + N)) ** 0.5) ** 0.5
    z_mm_0 = -((N + 2) ** 2 - (2 * gamma ** 2) / (k ** 2) - ((2 * gamma) / (k ** 2)) * (
                gamma ** 2 - 4 * k * (-1) ** n * (1 + N)) ** 0.5) ** 0.5
    z_pm_0 = +((N + 2) ** 2 - (2 * gamma ** 2) / (k ** 2) - ((2 * gamma) / (k ** 2)) * (
                gamma ** 2 - 4 * k * (-1) ** n * (1 + N)) ** 0.5) ** 0.5
    z_mp_0 = -((N + 2) ** 2 - (2 * gamma ** 2) / (k ** 2) + ((2 * gamma) / (k ** 2)) * (
                gamma ** 2 - 4 * k * (-1) ** n * (1 + N)) ** 0.5) ** 0.5

    n = 1
    z_pp_1 = +((N + 2) ** 2 - (2 * gamma ** 2) / (k ** 2) + ((2 * gamma) / (k ** 2)) * (
                gamma ** 2 - 4 * k * (-1) ** n * (1 + N)) ** 0.5) ** 0.5
    z_mm_1 = -((N + 2) ** 2 - (2 * gamma ** 2) / (k ** 2) - ((2 * gamma) / (k ** 2)) * (
                gamma ** 2 - 4 * k * (-1) ** n * (1 + N)) ** 0.5) ** 0.5
    z_pm_1 = +((N + 2) ** 2 - (2 * gamma ** 2) / (k ** 2) - ((2 * gamma) / (k ** 2)) * (
                gamma ** 2 - 4 * k * (-1) ** n * (1 + N)) ** 0.5) ** 0.5
    z_mp_1 = -((N + 2) ** 2 - (2 * gamma ** 2) / (k ** 2) + ((2 * gamma) / (k ** 2)) * (
                gamma ** 2 - 4 * k * (-1) ** n * (1 + N)) ** 0.5) ** 0.5

    fig, ax = plt.subplots(1, 1, figsize=(8, 2.5))

    scattersize = 70

    ax.scatter(0, 0, zorder=5, color="k", edgecolor="k", s=scattersize, linewidths=2)
    ax.scatter(np.pi, 0,zorder=5, color="k", edgecolor="k", s=scattersize, linewidths=2)
    ax.scatter(-np.pi, 0, zorder=5, color="k", edgecolor="k", s=scattersize, linewidths=2)
    #ax.scatter(0, z_mp_0,zorder=5, color="k", edgecolor="k", s=scattersize, linewidths=2)
    #ax.scatter(0, z_pp_0, zorder=5, color="k", edgecolor="k", s=scattersize, linewidths=2)
    #ax.scatter(0, z_pm_0, zorder=5, color="w", edgecolor="k", s=scattersize, linewidths=2)
    #ax.scatter(0, z_mm_0, zorder=5, color="w", edgecolor="k", s=scattersize, linewidths=2)

    xticks = [-2 * np.pi, -np.pi, 0, np.pi, 2 * np.pi]
    xtick_labels = [r"$-2\pi$", r"$-\pi$", r"$0$", r"$\pi$", r"$2\pi$"]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xtick_labels)
    ax.tick_params(axis="both", direction="in", labelsize=23)

    ax.plot(arg_variable, z_variable[::], color="gray", lw=0.7)
    ax.plot(arg_variable, -z_variable[::], color="gray", lw=0.7)
    ax.plot(arg_variable + 2 * np.pi, z_variable[::], color="gray", lw=0.7)
    ax.plot(arg_variable + 2 * np.pi, -z_variable[::], color="gray", lw=0.7)
    ax.plot(arg_variable - 2 * np.pi, z_variable[::], color="gray", lw=0.7)
    ax.plot(arg_variable - 2 * np.pi, -z_variable[::], color="gray", lw=0.7)
    ax.set_xlim([-2 * np.pi, 2 * np.pi])
    ax.set_ylim([-N * 1.02, N * 1.02])
    ax.set_xlabel('$\phi$', fontsize=28)
    ax.set_ylabel('$z$', fontsize=28)

    plt.tight_layout()
    plt.savefig("test.png", dpi=300)