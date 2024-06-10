from back_process import *
from functions import *
from jacobians import *


def zeta_01(x, A, sigma, k):
    return A * (np.exp((- x ** 2) / (2 * sigma ** 2))) * np.cos(k * x)


def zeta_02(x, A, sigma, k):
    return -A * np.sqrt(1 / (2 * sigma ** 2)) * x * (np.exp((- x ** 2) / (2 * sigma ** 2))) * np.sin(k * x)

def gaussian(x, A, x_0, sigma):
    return A * np.exp(- (x - x_0) ** 2 / (2 * sigma ** 2))

if __name__ == '__main__':
    disc = "C:"
    nu = "0.018"
    gamma_str = "0.185"
    sigma_str = "6.000"
    alpha_str = "6.524"
    directory = disc + '/mnustes_science/simulation_data/FD/PT_dimer/alpha=' + alpha_str + '/beta=1.000/mu=0.100/nu=' + nu + '/sigma=' + sigma_str + '/gamma=' + gamma_str
    Z_r_00 = np.loadtxt(directory + '/phi00/field_real_0.txt', delimiter=',')
    Z_i_00 = np.loadtxt(directory + '/phi00//field_img_0.txt', delimiter=',')
    params_00 = np.loadtxt(directory + '/phi00/parameters.txt', delimiter=',')

    Z_r_01 = np.loadtxt(directory + '/phi01/field_real_0.txt', delimiter=',')
    Z_i_01 = np.loadtxt(directory + '/phi01/field_img_0.txt', delimiter=',')
    X = np.loadtxt(directory + '/phi01/X.txt', delimiter=',')
    params_01 = np.loadtxt(directory + '/phi01/parameters.txt', delimiter=',')

    dx = X[1] - X[0]
    # N_resample = 4
    # [Z_r_00, Z_i_00, Z_r_01, Z_i_01], x_grid = resample([Z_r_00, Z_i_00, Z_r_01, Z_i_01], x_grid, N_resample)

    distances = [30.3, 30.4] #np.arange(30, 30.4, 0.1)
    Nx = len(X)
    real_guess = []
    imag_guess = []
    for d in distances:
        print("############ d = " + str(d) + " ############")
        X_L = - d / 2
        X_R = + d / 2
        J_L = np.argmin(np.abs(X - X_L))
        J_R = np.argmin(np.abs(X - X_R))
        J_center = np.argmin(np.abs(X))
        Delta_J_L = J_center - J_L
        Delta_J_R = J_R - J_center
        PHI_L = np.append(Z_r_00[Delta_J_L:], np.zeros(Delta_J_L)) + 1j * np.append(Z_i_00[Delta_J_L:], np.zeros(Delta_J_L))
        PHI_R = np.append(np.zeros(Delta_J_R), Z_r_01[:-Delta_J_R]) + 1j * np.append(np.zeros(Delta_J_R), Z_i_01[:-Delta_J_R])
        PHI_L = 1j * PHI_L
        PHI_R = 1j * PHI_R
        PHI = PHI_L + PHI_R
        fit_R = np.real(PHI)
        fit_I = np.imag(PHI)
        plt.close()
        [alpha, beta, gamma_0, mu, nu, sigma, phi] = params_00
        gamma_real = -(gaussian(X, gamma_0, -d / 2, sigma) - gaussian(X, gamma_0, d / 2, sigma))
        gamma_img = 0 * gamma_0 * np.sin(phi) * np.exp(- X ** 2 / (2 * sigma ** 2))
        gamma = [gamma_real, gamma_img]
        #plt.plot(X, fit_R, label="r")
        #plt.plot(X, fit_I, label="i")
        #plt.plot(X, gamma_real, c="k")
        #plt.legend()
        #plt.show()
        parameters = [alpha, beta, gamma, mu, nu]
        D2 = sparse_DD_neumann(Nx, dx)
        operators = np.array([D2])
        [F, G] = equations_FD("PDNLS", [fit_R, fit_I], [0], X, [0], parameters, operators)
        J = jacobians_FD("PDNLS", [fit_R, fit_I], [0], X, [0], parameters, operators)
        b = -np.append(F, G)
        dX = np.linalg.solve(J, b)
        res = np.dot(np.append(fit_R, fit_I) - dX, np.append(fit_R, fit_I) - dX)
        i = 0
        while i < 1000:
            fit_R = fit_R + dX[:Nx]
            fit_I = fit_I + dX[Nx:]
            [F, G] = equations_FD("PDNLS", [fit_R, fit_I], [0], X, [0], parameters, operators)
            J = jacobians_FD("PDNLS", [fit_R, fit_I], [0], X, [0], parameters, operators)
            b = -np.append(F, G)
            dX = np.linalg.solve(J, b)
            res = np.dot(dX, dX)
            print(res)
            i = i + 1
        real_guess.append(fit_R + dX[:Nx])
        imag_guess.append(fit_I + dX[Nx:])
    real_guess = np.array(real_guess)
    imag_guess = np.array(imag_guess)

    fig, ((ax1a, ax2a, ax3a)) = plt.subplots(nrows=1, ncols=3, figsize=(3, 1))

    pc_01a = ax1a.pcolor(X, distances, real_guess, cmap=parula_map)
    ax1a.set_ylabel("$d$", fontsize=15)
    ax1a.tick_params(axis="x", direction="in", labeltop=True, labelbottom=False, top=True, bottom=True)
    ax1a.tick_params(axis="y", direction="in", left=True, right=True)
    ax1a.grid(alpha=0.2, color="k")

    pc_02a = ax2a.pcolor(X, distances, imag_guess, cmap=parula_map)
    ax2a.set_yticklabels([])
    ax2a.tick_params(axis="x", direction="in", labeltop=True, labelbottom=False, top=True, bottom=True)
    ax2a.tick_params(axis="y", direction="in", left=True, right=True)
    ax2a.grid(alpha=0.2, color="k")

    pc_03a = ax3a.pcolor(X, distances, np.abs(real_guess + 1j * imag_guess), cmap=parula_map)
    ax3a.set_yticklabels([])
    ax3a.tick_params(axis="x", direction="in", labeltop=True, labelbottom=False, top=True, bottom=True)
    ax3a.tick_params(axis="y", direction="in", left=True, right=True)
    ax3a.grid(alpha=0.2, color="k")

    plt.show()
    plt.close()

