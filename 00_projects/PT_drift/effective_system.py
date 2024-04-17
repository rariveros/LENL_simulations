import matplotlib.pyplot as plt
import numpy as np

from functions import *
from back_process import *
from time_integrators import *
from scipy import integrate


def gaussian(x, A, x_0, sigma):
    return A * np.exp(- (x - x_0) ** 2 / (2 * sigma ** 2))


if __name__ == '__main__':
    disco = 'D:/'
    initial_dir_data = str(disco) + 'Users/mnustes_science/PT_fluids/mnustes_science/simulation_data'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    save_directory = directory + "/analysis"

    Z_r_00 = np.loadtxt(directory + '/phi00/field_real_0.txt', delimiter=',')
    Z_i_00 = np.loadtxt(directory + '/phi00//field_img_0.txt', delimiter=',')
    params_00 = np.loadtxt(directory + '/phi00/parameters.txt', delimiter=',')

    Z_r_01 = np.loadtxt(directory + '/phi01/field_real_0.txt', delimiter=',')
    Z_i_01 = np.loadtxt(directory + '/phi01/field_img_0.txt', delimiter=',')
    x_grid = np.loadtxt(directory + '/phi01/X.txt', delimiter=',')
    params_01 = np.loadtxt(directory + '/phi01/parameters.txt', delimiter=',')

    Nx = len(x_grid)
    distances = np.arange(2, 200, 1)
    PAR = []
    COF = []
    d = 150
    X_L = - d / 2
    X_R = + d / 2
    J_L = np.argmin(np.abs(x_grid - X_L))
    J_R = np.argmin(np.abs(x_grid - X_R))
    J_center = np.argmin(np.abs(x_grid))
    Delta_J_L = J_center - J_L
    Delta_J_R = J_R - J_center
    PHI_L = np.append(Z_r_00[Delta_J_L:], np.zeros(Delta_J_L)) + 1j * np.append(Z_i_00[Delta_J_L:], np.zeros(Delta_J_L))
    PHI_R = np.append(np.zeros(Delta_J_R), Z_r_01[:-Delta_J_R]) + 1j * np.append(np.zeros(Delta_J_R), Z_i_01[:-Delta_J_R])
    phi_02 = PHI_L + 1j * PHI_R
    phi_01 = PHI_L - 1j * PHI_R

    [alpha, beta, gamma_0, mu, nu, sigma, phi] = params_00
    gamma = gaussian(x_grid, gamma_0, -d / 2, sigma) - gaussian(x_grid, gamma_0, d / 2, sigma)
    plt.plot(x_grid, np.real(phi_01), color="purple", label="$\phi_{+}^{R}$")
    plt.plot(x_grid, np.imag(phi_01), color="green", label="$\phi_{+}^{I}$")
    plt.plot(x_grid, np.real(phi_02), color="purple", linestyle="--", label="$\phi_{-}^{R}$")
    plt.plot(x_grid, np.imag(phi_02), color="green", linestyle="--", label="$\phi_{-}^{I}$")

    #plt.plot(x_grid, np.real(PHI_R), color="b", label="$\phi_{R}^{R}$")
    #plt.plot(x_grid, np.imag(PHI_R), color="r", label="$\phi_{R}^{I}$")
    #plt.plot(x_grid, np.real(PHI_L), color="b", linestyle="--", label="$\phi_{L}^{R}$")
    #plt.plot(x_grid, np.imag(PHI_L), color="r", linestyle="--", label="$\phi_{L}^{I}$")
    plt.plot(x_grid, gamma, color="k", label="$\gamma(x)$")

    plt.legend(fontsize=15)
    plt.show()
    plt.close()
    for d in distances:
        X_L = - d / 2
        X_R = + d / 2
        J_L = np.argmin(np.abs(x_grid - X_L))
        J_R = np.argmin(np.abs(x_grid - X_R))
        J_center = np.argmin(np.abs(x_grid))
        Delta_J_L = J_center - J_L
        Delta_J_R = J_R - J_center

        PHI_L = np.append(Z_r_00[Delta_J_L:], np.zeros(Delta_J_L)) + 1j * np.append(Z_i_00[Delta_J_L:], np.zeros(Delta_J_L))
        PHI_R = np.append(np.zeros(Delta_J_R), Z_r_01[:-Delta_J_R]) + 1j * np.append(np.zeros(Delta_J_R), Z_i_01[:-Delta_J_R])

        [alpha, beta, gamma_0, mu, nu, sigma, phi] = params_00
        #print([alpha, beta, gamma_0, mu, nu, sigma, phi])
        sigma_w = np.sqrt(np.sqrt(nu * alpha) * (sigma / mu))


        phi_01 = PHI_L + 1j * PHI_R
        phi_02 = PHI_L - 1j * PHI_R
        gamma = gaussian(x_grid, gamma_0, -d / 2, sigma) - gaussian(x_grid, gamma_0, d / 2, sigma)

        Gamma_0 = integrate.trapezoid(phi_01 * phi_02, x_grid)
        Gamma_1 = integrate.trapezoid(phi_01 ** 2, x_grid)
        Gamma_2 = integrate.trapezoid(phi_02 ** 2, x_grid)

        ddx_phi_01 = np.append(0, np.append(0, np.diff(np.diff(phi_01))))
        ddx_phi_02 = np.append(0, np.append(0, np.diff(np.diff(phi_02))))

        Sigma_11 = integrate.simpson(phi_01 * ddx_phi_01, x_grid) / Gamma_1
        Sigma_12 = integrate.simpson(phi_01 * ddx_phi_02, x_grid) / Gamma_1
        Sigma_21 = integrate.simpson(phi_02 * ddx_phi_01, x_grid) / Gamma_2
        Sigma_22 = integrate.simpson(phi_02 * ddx_phi_02, x_grid) / Gamma_2

        Pi_11 = integrate.simpson(phi_01 * np.conjugate(phi_01) * gamma, x_grid) / Gamma_1
        Pi_12 = integrate.simpson(phi_01 * np.conjugate(phi_02) * gamma, x_grid) / Gamma_1
        Pi_21 = integrate.simpson(phi_02 * np.conjugate(phi_01) * gamma, x_grid) / Gamma_2
        Pi_22 = integrate.simpson(phi_02 * np.conjugate(phi_02) * gamma, x_grid) / Gamma_2

        a_L1 = (alpha / 2) * (Sigma_11 + Sigma_12 + Sigma_21 + Sigma_22)
        a_L2 = (alpha / 2) * (Sigma_11 - Sigma_12 + Sigma_21 - Sigma_22)
        b_L1 = (1 / 2) * (Pi_11 + Pi_12 + Pi_21 + Pi_22)
        b_L2 = (1 / 2) * (Pi_11 - Pi_12 + Pi_21 - Pi_22)
        a_R1 = (alpha / 2) * (-Sigma_11 - Sigma_12 + Sigma_21 + Sigma_22)
        a_R2 = (alpha / 2) * (-Sigma_11 + Sigma_12 + Sigma_21 - Sigma_22)
        b_R1 = (1 / 2) * (Pi_11 + Pi_12 - Pi_21 - Pi_22)
        b_R2 = (1 / 2) * (-Pi_11 + Pi_12 + Pi_21 - Pi_22)

        dJ = np.array([[- mu + np.imag(a_L1) + np.real(b_L1), -nu - np.real(a_L1) + np.imag(b_L1) , -np.real(a_R1) - np.imag(b_R1), -np.imag(a_R1) + np.real(b_R1)],
                       [+nu + np.real(a_L1) + np.imag(b_L1) , - mu + np.imag(a_L1) - np.real(b_L1), +np.imag(a_R1) + np.real(b_R1), -np.real(a_R1) + np.imag(b_R1)],
                       [-np.real(a_L2) - np.imag(b_L2)      ,  -np.imag(a_L2) + np.real(b_L2)     , - mu + np.imag(a_R2) + np.real(b_R2), -nu - np.real(a_R2) + np.imag(b_R2)],
                       [+np.imag(a_L2) + np.real(b_L2)      , -np.real(a_L2) + np.imag(b_L2)      , +nu + np.real(a_R2) + np.imag(b_R2) , - mu + np.imag(a_R2) - np.real(b_R2)]])
        eigenvalues, eigenvectors = np.linalg.eig(dJ)
        eig_R = np.real(eigenvalues)
        eig_I = np.imag(eigenvalues)
        #print(eig_I)
        for i in [0, 1, 2, 3]:
            if i == 0 and d == distances[0]:
                plt.scatter(d, eig_R[i], c="b", s=10, label="$\lambda_R$")
                plt.scatter(d, eig_I[i], c="r", s=10, label="$\lambda_I$")
            else:
                plt.scatter(d, eig_R[i], c="b", s=10)
                plt.scatter(d, eig_I[i], c="r", s=10)
        cof = [Gamma_1, Gamma_2, Pi_11, Pi_12, Pi_21, Pi_22, Sigma_11, Sigma_12, Sigma_21, Sigma_22, Gamma_0]
        parameters = [a_L1, a_L2, a_R1, a_R2, b_L1, b_L2, b_R1, b_R2]
        PAR.append(parameters)
        COF.append(cof)
    PAR = np.array(PAR)
    COF = np.array(COF)
    # plt.legend()
    plt.grid(alpha=0.2)
    plt.xticks(fontsize=15)
    plt.xlabel("$d$", fontsize=20)
    # plt.xlim(20, 40)
    plt.yticks(fontsize=15)
    plt.ylabel("$\lambda_R, \lambda_I$", fontsize=20)
    # plt.ylim(-30, 30)
    plt.legend(fontsize=15)
    plt.tight_layout()
    plt.savefig("eigen.png", dpi=200)
    plt.close()

    fig, axs = plt.subplots(nrows=3, ncols=2)
    axs[0, 0].plot(distances, np.real(COF[:, 0]), label="$\Gamma_{+}$")
    axs[0, 0].plot(distances, np.real(COF[:, 1]), label="$\Gamma_{-}$")
    axs[0, 0].plot(distances, np.real(COF[:, -1]), label="$\Gamma_{0}$")
    axs[0, 0].legend()
    axs[1, 0].plot(distances, np.imag(COF[:, 0]), label="$\Gamma_{+}$")
    axs[1, 0].plot(distances, np.imag(COF[:, 1]), label="$\Gamma_{-}$")
    axs[1, 0].legend()
    axs[0, 1].plot(distances, np.real(COF[:, 2]), label="$\Sigma_{++}$")
    axs[0, 1].plot(distances, np.real(COF[:, 3]), label="$\Sigma_{+-}$")
    axs[0, 1].plot(distances, np.real(COF[:, 4]), label="$\Sigma_{-+}$")
    axs[0, 1].plot(distances, np.real(COF[:, 5]), label="$\Sigma_{--}$")
    axs[0, 1].legend()
    axs[1, 1].plot(distances, np.imag(COF[:, 2]), label="$\Sigma_{++}$")
    axs[1, 1].plot(distances, np.imag(COF[:, 3]), label="$\Sigma_{+-}$")
    axs[1, 1].plot(distances, np.imag(COF[:, 4]), label="$\Sigma_{-+}$")
    axs[1, 1].plot(distances, np.imag(COF[:, 5]), label="$\Sigma_{--}$")
    axs[1, 1].legend()
    axs[2, 0].plot(distances, np.real(COF[:, 6]), label="$\Pi_{++}$")
    axs[2, 0].plot(distances, np.real(COF[:, 7]), label="$\Pi_{+-}$")
    axs[2, 0].plot(distances, np.real(COF[:, 8]), label="$\Pi_{-+}$")
    axs[2, 0].plot(distances, np.real(COF[:, 9]), label="$\Pi_{--}$")
    axs[2, 0].legend()
    axs[2, 1].plot(distances, np.imag(COF[:, 6]), label="$\Pi_{++}$")
    axs[2, 1].plot(distances, np.imag(COF[:, 7]), label="$\Pi_{+-}$")
    axs[2, 1].plot(distances, np.imag(COF[:, 8]), label="$\Pi_{-+}$")
    axs[2, 1].plot(distances, np.imag(COF[:, 9]), label="$\Pi_{--}$")
    axs[2, 1].legend()
    plt.show()
    plt.close()

    fig, axs = plt.subplots(nrows=2, ncols=2)
    axs[0, 0].plot(distances, np.real(PAR[:, 0]), label="$a_{L}^{1}$")
    axs[0, 0].plot(distances, np.real(PAR[:, 1]), label="$a_{L}^{2}$")
    axs[0, 0].plot(distances, np.real(PAR[:, 2]), label="$a_{R}^{1}$")
    axs[0, 0].plot(distances, np.real(PAR[:, 3]), label="$a_{R}^{2}$")
    axs[0, 0].legend()
    axs[1, 0].plot(distances, np.real(PAR[:, 4]), label="$b_{L}^{1}$")
    axs[1, 0].plot(distances, np.real(PAR[:, 5]), label="$b_{L}^{2}$")
    axs[1, 0].plot(distances, np.real(PAR[:, 6]), label="$b_{R}^{1}$")
    axs[1, 0].plot(distances, np.real(PAR[:, 7]), label="$b_{R}^{2}$")
    axs[1, 0].legend()
    axs[0, 1].plot(distances, np.imag(PAR[:, 0]), label="$a_{L}^{1}$")
    axs[0, 1].plot(distances, np.imag(PAR[:, 1]), label="$a_{L}^{2}$")
    axs[0, 1].plot(distances, np.imag(PAR[:, 2]), label="$a_{R}^{1}$")
    axs[0, 1].plot(distances, np.imag(PAR[:, 3]), label="$a_{R}^{2}$")
    axs[0, 1].legend()
    axs[1, 1].plot(distances, np.imag(PAR[:, 4]), label="$b_{L}^{1}$")
    axs[1, 1].plot(distances, np.imag(PAR[:, 5]), label="$b_{L}^{2}$")
    axs[1, 1].plot(distances, np.imag(PAR[:, 6]), label="$b_{R}^{1}$")
    axs[1, 1].plot(distances, np.imag(PAR[:, 7]), label="$b_{R}^{2}$")
    axs[1, 1].legend()
    plt.tight_layout()
    plt.savefig("coeficientes.png", dpi=300)
    plt.close()






