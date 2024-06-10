from functions import *
from back_process import *
from time_integrators import *


if __name__ == '__main__':
    x_grid = np.arange(-150, 150, 0.1)
    dx = x_grid[1] - x_grid[0]
    alpha = 26.096
    beta = 1
    sigma_i = 15
    nu = -0.13
    mu = 0.075
    gamma_0 = 0.2

    X = np.arange(-12, 12, 0.25)
    F = []
    for i in range(len(X)):
        gammaR_i = gamma_0 * (np.exp(-X[i] ** 2 / (2 * sigma_i ** 2)))
        gammaI_i = gamma_0 * (np.exp(-X[i] ** 2 / (2 * sigma_i ** 2)))
        phaseR_i = 0.5 * np.arccos(mu / gammaR_i)
        phaseI_i = 0.5 * np.arccos(mu / gammaI_i)
        deltaR_i = 1 * (- nu + np.sqrt(gammaR_i ** 2 - mu ** 2))
        deltaI_i = 1 * (- nu + np.sqrt(gammaI_i ** 2 - mu ** 2))
        AR_i = np.sqrt(2 * deltaR_i) * np.cos(phaseR_i)
        AI_i = np.sqrt(2 * deltaI_i) * np.sin(phaseI_i)
        phi_01_i = (AR_i / np.cosh(np.sqrt(deltaR_i / alpha) * (x_grid - (X[i]))))
        phi_02_i = - (AI_i / np.cosh(np.sqrt(deltaI_i / alpha) * (x_grid - (X[i] + 5))))
        #plt.plot(x_grid, phi_01_i, color="b")
        #plt.plot(x_grid, phi_02_i, color="r")
        #plt.plot(x_grid, (phi_01_i ** 2 + phi_02_i ** 2) ** 0.5, color="k")
        #plt.plot(x_grid, gamma_0 * (np.exp(-x_grid ** 2 / (2 * sigma_i ** 2))), color="k", linestyle="--")
        #plt.hlines(0.075, -120, 120, colors="k", linestyles="--")
        #plt.hlines(0., -120, 120, colors="k")
        #plt.grid(alpha=0.2)
        #plt.show()
        #plt.close()
        Dphi_01_i = np.append(np.diff(phi_01_i) / (dx), 0)
        Dphi_02_i = np.append(np.diff(phi_02_i) / (dx), 0)
        f_i = (-(alpha / 2) * (Dphi_01_i ** 2 + Dphi_02_i ** 2) + (beta / 4) * (phi_01_i ** 2 + phi_02_i ** 2) ** 2 + (nu / 2) * (phi_01_i ** 2 + phi_02_i ** 2) + gamma_0 * (np.exp(-x_grid ** 2/(2 * sigma_i ** 2))) * phi_01_i * phi_02_i) / (2 * mu)
        F_i = integrate.simpson(f_i, x_grid)
        F.append(F_i)
    F = np.array(F)
    dF = np.diff(F) / (X[1] - X[0])
    ddF = np.diff(dF) / (X[1] - X[0])
    argmin = np.argmax(np.real(F))
    print(X[argmin])
    plt.plot(X, np.real(F))
    #plt.plot(X, np.imag(ddF))
    plt.show()
    #X, Y = np.meshgrid(b, A)
    #fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    #surf = ax.plot_surface(X, Y, F, cmap=parula_map, linewidth=0, antialiased=False)
    #ax.set_xlabel("$b$")
    #ax.set_ylabel("$A$")
    #plt.show()