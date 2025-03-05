from back_process import *
from functions import *
from jacobians import *

if __name__ == '__main__':
    directory = 'D:/mnustes_science/simulation_data/FD/soliton_control/alpha=6.524/beta=1.000/mu=0.075/nu=-0.150/sigma=15.000/gamma=0.180'

    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')[-1, :]
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')[-1, :]
    X = np.loadtxt(directory + '/X.txt', delimiter=',')
    [alpha, beta, gamma_0, mu, nu] = np.loadtxt(directory + '/parameters.txt', delimiter=',')
    sigma = 15
    phi = 0
    Nx = len(X)
    dx = X[1] - X[0]
    D2 = sparse_DD_neumann(Nx, dx)
    operators = np.array([D2])

    nu = -0.151
    gamma_real = gamma_0 * np.exp(- X ** 2 / (2 * sigma ** 2))
    gamma_img = gamma_0 * np.exp(- X ** 2 / (2 * sigma ** 2)) * 0
    gamma = [gamma_real, gamma_img]
    parameters = [alpha, beta, gamma, mu, nu]

    # Initialize state vector
    U = np.concatenate([Z_r, Z_i])  # Vector unificado
    Nx = len(X)
    tol = 1e-6  # Tolerancia más estricta para estabilidad
    max_iter = 100  # Aumentamos iteraciones
    damping_factor = 1.0  # Inicialmente 1 (Newton estándar)

    # Primera evaluación del sistema
    F, G = equations_FD("PDNLS", [U[:Nx], U[Nx:]], [0], X, [0], parameters, operators)
    J = jacobians_FD("PDNLS", [U[:Nx], U[Nx:]], [0], X, [0], parameters, operators)

    # Regularización del Jacobiano para evitar singularidades
    J += np.eye(J.shape[0]) * 1e-6  # Pequeño término diagonal

    # Resolver sistema J * dX = -F
    b = -np.concatenate([F, G])

    try:
        dX = np.linalg.solve(J, b)  # Más estable que inv(J) @ b
    except np.linalg.LinAlgError:
        print("Warning: Jacobian singular! Using least squares solution.")
        dX = -np.linalg.lstsq(J, b, rcond=None)[0]  # Usa mínimos cuadrados si J es singular

    # Calcular norma del paso
    res = np.linalg.norm(dX)
    print(f"Initial Residual: {res:.2e}")

    # Iteración de Newton
    i = 0
    while res > tol and i < max_iter:
        # Aplicar damping si es necesario
        U_new = U + damping_factor * dX

        # Evaluar ecuaciones con nuevo U
        F, G = equations_FD("PDNLS", [U_new[:Nx], U_new[Nx:]], [0], X, [0], parameters, operators)
        J = jacobians_FD("PDNLS", [U_new[:Nx], U_new[Nx:]], [0], X, [0], parameters, operators)

        # Regularización
        J += np.eye(J.shape[0]) * 1e-6

        b = -np.concatenate([F, G])

        try:
            dX_new = np.linalg.solve(J, b)
        except np.linalg.LinAlgError:
            dX_new = -np.linalg.lstsq(J, b, rcond=None)[0]

            # Nuevo residuo
        res_new = np.linalg.norm(dX_new)

        # Damping adaptativo: reducir el paso si el residuo crece
        if res_new > res:
            damping_factor *= 0.5
            print(f"Reducing step size: λ = {damping_factor:.2f}")
        else:
            U = U_new
            dX = dX_new
            res = res_new

        print(f"Iteration {i + 1}: Residual = {res:.2e}")
        i += 1

    # Extraer partes real e imaginaria
    fit_R, fit_I = U[:Nx], U[Nx:]

    # Gráficos de comparación
    plt.plot(X, Z_r, color="b", label="Real (True)")
    plt.plot(X, Z_i, color="r", label="Imag (True)")
    plt.plot(X, fit_R, '--', color="b", label="Real (Optimized)")
    plt.plot(X, fit_I, '--', color="r", label="Imag (Optimized)")
    plt.legend()
    plt.show()