from functions import *
from back_process import *
from time_integrators import *
from scipy.signal import argrelextrema

if __name__ == '__main__':
    x_grid = np.arange(-75, 75, 0.1)
    dx = x_grid[1] - x_grid[0]
    alpha = 26.096
    beta = 1
    sigma_i = 17.5
    nu = -0.1
    mu = 0.075
    gammas = np.arange(0.1, 0.25, 0.01)

    X = np.arange(-10, 10, 0.1)
    Y = np.arange(-10, 10, 0.1)
    Xs = []
    Ys = []
    nus_scat = []
    for gamma_0 in gammas:
        print(gamma_0)
        F = []
        for i in range(len(X)):
            Fi = []
            for j in range(len(Y)):
                gammaR_ij = gamma_0 * (np.exp(-X[i] ** 2 / (2 * sigma_i ** 2)))
                gammaI_ij = gamma_0 * (np.exp(-Y[j] ** 2 / (2 * sigma_i ** 2)))
                phaseR_ij = 0.5 * np.arccos(mu / gammaR_ij)
                phaseI_ij = 0.5 * np.arccos(mu / gammaI_ij)
                deltaR_ij = 1 * (- nu + np.sqrt(gammaR_ij ** 2 - mu ** 2))
                deltaI_ij = 1 * (- nu + np.sqrt(gammaI_ij ** 2 - mu ** 2))
                AR_i = np.sqrt(2 * deltaR_ij) * np.cos(phaseR_ij)
                AI_i = - np.sqrt(2 * deltaI_ij) * np.sin(phaseI_ij)
                phi_01_i = (AR_i / np.cosh(np.sqrt(deltaR_ij / alpha) * (x_grid - X[i])))
                phi_02_i = (AI_i / np.cosh(np.sqrt(deltaI_ij / alpha) * (x_grid - Y[j])))
                Dphi_01_i = np.append(np.diff(phi_01_i) / (dx), 0)
                Dphi_02_i = np.append(np.diff(phi_02_i) / (dx), 0)
                f_ij = (-(alpha / 2) * (Dphi_01_i ** 2 + Dphi_02_i ** 2) + (beta / 4) * (phi_01_i ** 2 + phi_02_i ** 2) ** 2 + (nu / 2) * (phi_01_i ** 2 + phi_02_i ** 2) + gamma_0 * (np.exp(-x_grid ** 2/(2 * sigma_i ** 2))) * phi_01_i * phi_02_i) / (2 * mu)
                F_ij = integrate.simpson(f_ij, x_grid)
                Fi.append(F_ij)
            Fi = np.array(Fi)
            F.append(Fi)
        F = np.array(F)
        mins_i = []
        for i in range(len(X)):
            for j in range(len(Y)):
                if argrelextrema(F[i, :], np.less)[0][0] == j and argrelextrema(F[:, j], np.less)[0][0] == i:
                    Xs.append(X[i])
                    Ys.append(-Y[j])
                    nus_scat.append(gamma_0)

    plt.scatter(nus_scat, Xs, color="r")
    plt.scatter(nus_scat, Ys, color="b")
    plt.show()