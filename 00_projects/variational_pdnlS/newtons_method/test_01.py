from back_process import *
from functions import *
from jacobians import *


def zeta_01(x, A, sigma, k):
    return A * (np.exp((- x ** 2) / (2 * sigma ** 2))) * np.cos(k * x)


def zeta_02(x, A, sigma, k):
    return -A * np.sqrt(1 / (2 * sigma ** 2)) * x * (np.exp((- x ** 2) / (2 * sigma ** 2))) * np.sin(k * x)

if __name__ == '__main__':
    directory = 'C:/mnustes_science/simulation_data/FD/variational_pdnlS/localized_patterns/alpha=6.524/beta=1.000/mu=0.100/nu=0.018/sigma=6.000/gamma=0.185'

    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    X = np.loadtxt(directory + '/X.txt', delimiter=',')
    [alpha, beta, gamma_0, mu, nu, sigma, phi] = np.loadtxt(directory + '/parameters.txt', delimiter=',')
    Nx = len(X)
    dx = X[1] - X[0]
    D2 = sparse_DD_neumann(Nx, dx)
    operators = np.array([D2])

    param_01, cov_01 = curve_fit(zeta_01, X, Z_r, bounds=[(0.1, 0, 0.01), (0.5, 25, 0.2)])
    #A_01 = param_01[0]
    #sigma_01 = param_01[1]
    #k_01 = param_01[2]

    A_01 = 0.5
    sigma_01 = sigma
    k_01 = np.sqrt(nu / alpha)

    fit_R = zeta_01(X, A_01, sigma_01, k_01)
    fit_I = zeta_02(X, A_01, sigma_01, k_01)
    fit_R0 = fit_R
    fit_I0 = fit_I
    plt.plot(X, fit_R0, linestyle='dotted', color="b")
    plt.plot(X, fit_I0, linestyle='dotted', color='r')

    gamma_real = gamma_0 * np.cos(phi) * np.exp(- X ** 2 / (2 * sigma ** 2))
    gamma_img = gamma_0 * np.sin(phi) * np.exp(- X ** 2 / (2 * sigma ** 2))
    gamma = [gamma_real, gamma_img]
    parameters = [alpha, beta, gamma, mu, nu]


    [F, G] = equations_FD("PDNLS", [fit_R, fit_I], [0], X, [0], parameters, operators)
    J = jacobians_FD("PDNLS", [fit_R, fit_I], [0], X, [0], parameters, operators)
    b = -np.append(F, G)
    dX = np.linalg.solve(J, b)
    res = np.dot(np.append(fit_R, fit_I) - dX, np.append(fit_R, fit_I) - dX)
    print(res)
    i = 0
    while res < 100 and i < 10:
        fit_R = fit_R + dX[:Nx]
        fit_I = fit_I + dX[Nx:]
        [F, G] = equations_FD("PDNLS", [fit_R, fit_I], [0], X, [0], parameters, operators)
        J = jacobians_FD("PDNLS", [fit_R, fit_I], [0], X, [0], parameters, operators)
        b = -np.append(F, G)
        dX = np.linalg.solve(J, b)
        res = np.dot(dX, dX)
        print(res)
        i = i + 1
    plt.plot(X, Z_r, color="b")
    plt.plot(X, Z_i, color='r')
    plt.plot(X, fit_R + dX[:Nx], '--', color="b")
    plt.plot(X, fit_I + dX[Nx:], '--', color='r')
    plt.show()

