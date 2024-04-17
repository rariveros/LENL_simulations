import matplotlib.pyplot as plt
import numpy as np

from functions import *
from back_process import *
from time_integrators import *
from scipy import integrate

def gaussian(x, A, x_0, sigma):
    return A * np.exp(-(x - x_0) ** 2 / (2 * sigma ** 2))

if __name__ == '__main__':
    x_grid = np.arange(-300, 300, 0.02)
    A = 0.7
    alpha = 6.5
    sigma = 16 * np.sqrt(alpha)
    distances = np.arange(1, 150, 1)
    nu = 0.1
    mu = 0.1
    gamma_0 = 0.24
    sigma_w = np.sqrt(np.sqrt(nu*alpha) * sigma / mu)

    PAR = []
    for d in distances:
        phi_01 = gaussian(x_grid, A, -d/2, sigma_w) * np.cos(-np.sqrt(nu / alpha) * (x_grid + d/2))
        phi_02 = gaussian(x_grid, A, +d/2, sigma_w) * np.cos(np.sqrt(nu / alpha) * (x_grid - d/2))
        if d == distances[0]:
            plt.plot(x_grid, phi_01)
            plt.plot(x_grid, phi_02)
            plt.show()
            plt.close()
        gamma = -(gaussian(x_grid, gamma_0, -d/2, sigma) - gaussian(x_grid, gamma_0, d/2, sigma))

        Gamma_1 = integrate.simpson(phi_01 ** 2, x_grid)
        Gamma_2 = integrate.simpson(phi_02 ** 2, x_grid)

        ddx_phi_01 = np.append(np.append(0, np.diff(np.diff(phi_01))), 0)
        ddx_phi_02 = np.append(np.append(0, np.diff(np.diff(phi_02))), 0)

        Sigma_11 = integrate.simpson(phi_01 * ddx_phi_01, x_grid)
        Sigma_12 = integrate.simpson(phi_01 * ddx_phi_02, x_grid)
        Sigma_21 = integrate.simpson(phi_02 * ddx_phi_01, x_grid)
        Sigma_22 = integrate.simpson(phi_02 * ddx_phi_02, x_grid)

        Delta_1111 = integrate.simpson(phi_01 ** 4 * phi_02 ** 0, x_grid)
        Delta_1112 = integrate.simpson(phi_01 ** 3 * phi_02 ** 1, x_grid)
        Delta_1122 = integrate.simpson(phi_01 ** 2 * phi_02 ** 2, x_grid)
        Delta_1222 = integrate.simpson(phi_01 ** 1 * phi_02 ** 3, x_grid)
        Delta_2222 = integrate.simpson(phi_01 ** 0 * phi_02 ** 4, x_grid)

        Eta_1 = integrate.simpson(phi_01 * phi_02 * gamma, x_grid)
        Eta_2 = integrate.simpson(phi_01 * phi_01 * gamma, x_grid)
        Eta_3 = integrate.simpson(phi_02 * phi_02 * gamma, x_grid)

        dJ = np.array([[-mu + Eta_2 / Gamma_1, nu + Sigma_11 / Gamma_1, + Eta_1 / Gamma_1, + Sigma_12 / Gamma_1],
                       [- (nu + Sigma_11 / Gamma_1), -mu - Eta_2 / Gamma_1, - Sigma_12 / Gamma_1, - Eta_1 / Gamma_1],
                       [+Eta_1 / Gamma_2, +Sigma_21 / Gamma_2, -mu + Eta_3 / Gamma_2, nu + Sigma_22 / Gamma_2],
                       [- Sigma_21 / Gamma_2, - Eta_1 / Gamma_2, - (nu + Sigma_22 / Gamma_2), -mu - Eta_3 / Gamma_2]])
        eigenvalues, eigenvectors = np.linalg.eig(dJ)
        eig_R = np.real(eigenvalues)
        eig_I = np.imag(eigenvalues)
        for i in range(len(eig_R)):
            if i == 1 and d==distances[0]:
                plt.scatter(d, eig_R[i], c="b", s=10, label="$\lambda_R$")
                plt.scatter(d, eig_I[i], c="r", s=10, label="$\lambda_L$")
            else:
                plt.scatter(d, eig_R[i], c="b", s=10)
                plt.scatter(d, eig_I[i], c="r", s=10)

        parameters = [Gamma_1, Gamma_2, Sigma_11, Sigma_12, Sigma_21, Sigma_22, Delta_1111, Delta_1112, Delta_1122, Delta_1222, Delta_2222, Eta_1]
        PAR.append(parameters)
    PAR = np.array(PAR)
    plt.legend(fontsize=15)
    plt.grid(alpha=0.2)
    plt.xlabel("$a$", fontsize=20)
    plt.ylabel("$\lambda_R, \lambda_I$", fontsize=20)
    plt.show()
    plt.savefig("eigen.png", dpi=200)
    plt.close()

    fig, axs = plt.subplots(nrows=2, ncols=2)
    axs[0, 0].plot(distances, PAR[:, 0], label="$\Gamma_{+}$")
    axs[0, 0].plot(distances, PAR[:, 1], label="$\Gamma_{-}$")
    axs[0, 0].legend()
    axs[1, 0].plot(distances, PAR[:, 2], label="$\Sigma_{++}$")
    axs[1, 0].plot(distances, PAR[:, 3], label="$\Sigma_{+-}$")
    axs[1, 0].plot(distances, PAR[:, 4], label="$\Sigma_{-+}$")
    axs[1, 0].plot(distances, PAR[:, 5], label="$\Sigma_{--}$")
    axs[1, 0].legend()
    axs[0, 1].plot(distances, PAR[:, 6], label="$\Delta_{1}$")
    axs[0, 1].plot(distances, PAR[:, 7], label="$\Delta_{2}$")
    axs[0, 1].plot(distances, PAR[:, 8], label="$\Delta_{3}$")
    axs[0, 1].plot(distances, PAR[:, 9], label="$\Delta_{4}$")
    axs[0, 1].plot(distances, PAR[:, 10], label="$\Delta_{5}$")
    axs[0, 1].legend()
    axs[1, 1].plot(distances, PAR[:, 11], label="$\Pi_{+-}$")
    axs[1, 1].legend()
    plt.savefig("coeficientes.png", dpi=300)