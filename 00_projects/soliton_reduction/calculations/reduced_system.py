from functions import *
from back_process import *
from time_integrators import *

def phis(alpha, beta, nu, mu, gamma, sigma, X, Y, x_grid, dx):
    gamma_01 = gamma * np.exp(- X ** 2 / (2 * sigma ** 2))
    gamma_02 = gamma * np.exp(- Y ** 2 / (2 * sigma ** 2))
    delta_01 = - nu + np.sqrt(gamma_01 ** 2 - mu ** 2)
    delta_02 = - nu + np.sqrt(gamma_02 ** 2 - mu ** 2)
    theta_01 = 0.5 * np.arccos(mu / gamma_01)
    theta_02 = 0.5 * np.arccos(mu / gamma_02)
    phi_01 = (1 / (np.cosh(np.sqrt(delta_01 / alpha) * (x_grid - X)))) * np.sqrt(2 * delta_01) * np.cos(theta_01)
    phi_02 = - (1 / (np.cosh(np.sqrt(delta_02 / alpha) * (x_grid - Y)))) * np.sqrt(2 * delta_02) * np.sin(theta_02)
    Dphi_01 = np.append(np.diff(phi_01) / dx, 0)
    Dphi_02 = np.append(np.diff(phi_02) / dx, 0)
    return [phi_01, phi_02, Dphi_01, Dphi_02]


if __name__ == '__main__':
    [alpha, beta, nu, mu, gamma, sigma] = [1, 1, -0.2, 0.075, 0.18, 15]
    dx = 0.01
    x_grid = np.arange(-50, 50, dx)
    Xs = np.arange(-1, 1, 0.1)
    Ys = np.arange(-1, 1, 0.1)
    gamma_x = gamma * np.exp(- x_grid ** 2 / (2 * sigma ** 2))
    F_01 = []
    F_02 = []
    I1 = []
    for X in Xs:
        F_01_i = []
        F_02_i = []
        I1_i = []
        for Y in Ys:
            [phi_01, phi_02, Dphi_01, Dphi_02] = phis(alpha, beta, nu, mu, gamma, sigma, X, Y, x_grid, dx)
            i1 = 0.5 * phi_01 * phi_02
            i2a = 0.5 * (alpha * Dphi_01 ** 2 - nu * phi_01 ** 2 - (beta / 2) * phi_01 ** 4)
            i2b = 0.5 * (alpha * Dphi_02 ** 2 - nu * phi_02 ** 2 - (beta / 2) * phi_02 ** 4)
            i3 = (beta / 2) * phi_01 ** 2 * phi_02 ** 2 + gamma_x * phi_01 * phi_02
            F_01_ij = integrate.simpson(- 2 * mu * i1 + i2b + i3, x_grid)
            F_02_ij = integrate.simpson(- 2 * mu * i1 - i2b - i3, x_grid)
            I_01_ij = integrate.simpson(i1, x_grid)
            F_01_i.append(F_01_ij)
            F_02_i.append(F_02_ij)
            I1_i.append(I_01_ij)
        F_01.append(F_01_i)
        F_02.append(F_02_i)
        I1.append(I1_i)


