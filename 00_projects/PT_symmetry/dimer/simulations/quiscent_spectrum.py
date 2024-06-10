import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *
from numpy.fft import fft, ifft

def gaussian(x, A, x_0, sigma):
    return A * np.exp(- (x - x_0) ** 2 / (2 * sigma ** 2))


def resample(old_arrays, old_x, N_resample):
    new_arrays = []
    for i in range(len(old_arrays)):
        [new_array_i, new_x] = signal.resample(old_arrays[i], int(len(old_x) * (N_resample)), old_x)
        new_arrays.append(new_array_i)
    return new_arrays, new_x


if __name__ == '__main__':

    eq = 'PT_dimer'
    t_rate = 1
    dt = 1
    T = 5000

    disco = 'C:/'
    initial_dir_data = str(disco) + 'Users/mnustes_science/PT_fluids/mnustes_science/simulation_data'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    Z_r_00 = np.loadtxt(directory + '/phi00/field_real_0.txt', delimiter=',')
    Z_i_00 = np.loadtxt(directory + '/phi00//field_img_0.txt', delimiter=',')
    params_00 = np.loadtxt(directory + '/phi00/parameters.txt', delimiter=',')

    Z_r_01 = np.loadtxt(directory + '/phi01/field_real_0.txt', delimiter=',')
    Z_i_01 = np.loadtxt(directory + '/phi01/field_img_0.txt', delimiter=',')
    x_grid = np.loadtxt(directory + '/phi01/X.txt', delimiter=',')
    params_01 = np.loadtxt(directory + '/phi01/parameters.txt', delimiter=',')

    dx = x_grid[1] - x_grid[0]
    #N_resample = 4
    #[Z_r_00, Z_i_00, Z_r_01, Z_i_01], x_grid = resample([Z_r_00, Z_i_00, Z_r_01, Z_i_01], x_grid, N_resample)

    distances = np.arange(3, 60, 0.25)
    PHI_Rs = []
    PHI_Ls = []
    times = []
    UR_Rs = []
    UR_Is = []
    UL_Rs = []
    UL_Is = []
    for d in distances:
        print("############ d = " + str(d) + " ############")
        X_L = - d / 2
        X_R = + d / 2
        J_L = np.argmin(np.abs(x_grid - X_L))
        J_R = np.argmin(np.abs(x_grid - X_R))
        J_center = np.argmin(np.abs(x_grid))
        Delta_J_L = J_center - J_L
        Delta_J_R = J_R - J_center
        PHI_L = np.append(Z_r_00[Delta_J_L:], np.zeros(Delta_J_L)) + 1j * np.append(Z_i_00[Delta_J_L:], np.zeros(Delta_J_L))
        PHI_R = np.append(np.zeros(Delta_J_R), Z_r_01[:-Delta_J_R]) + 1j * np.append(np.zeros(Delta_J_R), Z_i_01[:-Delta_J_R])
        PHI_L = 1j * PHI_L
        PHI_R = 1j * PHI_R

        PHI = PHI_L + PHI_R
        phi_01 = PHI_L - 1j * PHI_R
        phi_02 = PHI_L + 1j * PHI_R

        if d == 320:
            plt.plot(x_grid, np.real(PHI), color="b")
            plt.plot(x_grid, np.imag(PHI), color="r")
            plt.plot(x_grid, np.abs(PHI), color="k")
            plt.hlines(0., -120, 120, colors="k", alpha=0.4)
            plt.grid(alpha=0.2)
            plt.show()
            plt.close()
        [alpha, beta, gamma_0, mu, nu, sigma, phi] = params_00
        gamma = -(gaussian(x_grid, gamma_0, -d / 2, sigma) - gaussian(x_grid, gamma_0, d / 2, sigma))

        Gamma_0 = integrate.simpson(phi_01 * phi_02, x_grid)
        Gamma_1 = integrate.simpson(phi_01 ** 2, x_grid)
        Gamma_2 = integrate.simpson(phi_02 ** 2, x_grid)

        ddx_phi_01 = alpha * np.append(0, np.append(np.diff(np.diff(phi_01)) / (dx ** 2), 0))
        ddx_phi_02 = alpha * np.append(0, np.append(np.diff(np.diff(phi_02)) / (dx ** 2), 0))

        Sigma_11 = integrate.simpson(phi_01 * ddx_phi_01, x_grid) / Gamma_1
        Sigma_12 = integrate.simpson(phi_01 * ddx_phi_02, x_grid) / Gamma_1
        Sigma_21 = integrate.simpson(phi_02 * ddx_phi_01, x_grid) / Gamma_2
        Sigma_22 = integrate.simpson(phi_02 * ddx_phi_02, x_grid) / Gamma_2

        Pi_11 = integrate.simpson(phi_01 * np.conjugate(phi_01) * gamma, x_grid) / Gamma_1
        Pi_12 = integrate.simpson(phi_01 * np.conjugate(phi_02) * gamma, x_grid) / Gamma_1
        Pi_21 = integrate.simpson(phi_02 * np.conjugate(phi_01) * gamma, x_grid) / Gamma_2
        Pi_22 = integrate.simpson(phi_02 * np.conjugate(phi_02) * gamma, x_grid) / Gamma_2

        Delta_11 = (1) * integrate.simpson(phi_01 * np.abs(phi_02) ** 2 * phi_02, x_grid) / Gamma_1
        Delta_21 = (2) * integrate.simpson(phi_01 * np.abs(phi_02) ** 2 * phi_01, x_grid) / Gamma_1
        Delta_31 = (2) * integrate.simpson(phi_01 * np.abs(phi_01) ** 2 * phi_02, x_grid) / Gamma_1
        Delta_41 = (1) * integrate.simpson(phi_01 * np.abs(phi_01) ** 2 * phi_01, x_grid) / Gamma_1
        Delta_51 = (1) * integrate.simpson(phi_01 * phi_02 ** 2 * np.conjugate(phi_01), x_grid) / Gamma_1
        Delta_61 = (1) * integrate.simpson(phi_01 * phi_01 ** 2 * np.conjugate(phi_02), x_grid) / Gamma_1

        Delta_12 = (1) * integrate.simpson(phi_02 * np.abs(phi_02) ** 2 * phi_02, x_grid) / Gamma_2
        Delta_22 = (2) * integrate.simpson(phi_02 * np.abs(phi_02) ** 2 * phi_01, x_grid) / Gamma_2
        Delta_32 = (2) * integrate.simpson(phi_02 * np.abs(phi_01) ** 2 * phi_02, x_grid) / Gamma_2
        Delta_42 = (1) * integrate.simpson(phi_02 * np.abs(phi_01) ** 2 * phi_01, x_grid) / Gamma_2
        Delta_52 = (1) * integrate.simpson(phi_02 * phi_02 ** 2 * np.conjugate(phi_01), x_grid) / Gamma_2
        Delta_62 = (1) * integrate.simpson(phi_02 * phi_01 ** 2 * np.conjugate(phi_02), x_grid) / Gamma_2
        #dJ = np.array([[- mu + np.imag(Sigma_11) + np.real(Pi_11), + nu + np.real(Sigma_11) + np.imag(Pi_11),
        #            + np.imag(Sigma_12) + np.real(Pi_12), + np.real(Sigma_12) + np.imag(Pi_12)],
        #           [- nu - np.real(Sigma_11) + np.imag(Pi_11), - mu + np.imag(Sigma_11) - np.real(Pi_11),
        #            - np.real(Sigma_12) + np.imag(Pi_12), np.imag(Sigma_12) - np.real(Pi_12)],
        #           [np.imag(Sigma_21) + np.real(Sigma_21), np.real(Sigma_21) + np.imag(Pi_21),
        #            - mu + np.imag(Sigma_22) + np.real(Pi_22), + nu + np.real(Sigma_22) + np.imag(Pi_22)],
        #           [-np.real(Sigma_21) + np.imag(Pi_21), np.imag(Sigma_21) - np.real(Pi_21),
        #            - nu - np.real(Sigma_22) + np.imag(Pi_22), - mu + np.imag(Sigma_22) - np.real(Pi_22)]])
        dJ = np.array([[- mu + np.real(Pi_11) + np.imag(Sigma_11), + nu + np.imag(Pi_11) + np.real(Sigma_11), + np.real(Pi_12) + np.imag(Sigma_12)     , + np.imag(Pi_12) + np.real(Sigma_12)],
                       [- nu + np.imag(Pi_11) - np.real(Sigma_11), - mu - np.real(Pi_11) + np.imag(Sigma_11), + np.imag(Pi_12) - np.real(Sigma_12)     , - np.real(Pi_12) + np.imag(Sigma_12)],
                       [np.real(Pi_21) + np.imag(Sigma_21)       , + np.imag(Pi_21) + np.real(Sigma_21)       , - mu + np.real(Pi_22) + np.imag(Sigma_22), + nu + np.imag(Pi_22) + np.real(Sigma_22)],
                       [np.imag(Pi_21) - np.real(Sigma_21)       , - np.real(Pi_21) + np.imag(Sigma_21)      , - nu + np.imag(Pi_22) - np.real(Sigma_22), - mu - np.real(Pi_22) + np.imag(Sigma_22)]])
        eigenvalues, eigenvectors = np.linalg.eig(dJ)
        eig_R = np.real(eigenvalues)
        eig_I = np.imag(eigenvalues)
        for i in [0, 1, 2, 3]:
            if i == 0 and d == distances[0]:
                plt.scatter(d, eig_R[i], c="b", s=10, label="$\lambda_R$")
                plt.scatter(d, eig_I[i], c="r", s=10, label="$\lambda_I$")
            else:
                plt.scatter(d, eig_R[i], c="b", s=10)
                plt.scatter(d, eig_I[i], c="r", s=10)

        sigma = [Sigma_11, Sigma_12, Sigma_21, Sigma_22]
        pi = [Pi_11, Pi_12, Pi_21, Pi_22]
        delta = [Delta_11, Delta_21, Delta_31, Delta_41, Delta_51, Delta_61, Delta_12, Delta_22, Delta_32, Delta_42, Delta_52, Delta_62]
        parameters = [alpha, beta, mu, nu, sigma, pi, delta]

    plt.grid(alpha=0.2)
    plt.xticks(fontsize=15)
    plt.xlabel("$d$", fontsize=20)
    plt.yticks(fontsize=15)
    plt.ylabel("$\lambda_R, \lambda_I$", fontsize=20)
    plt.hlines(0., 0, 100, colors="k", alpha=0.4)
    plt.xlim(distances[0], distances[-1])
    #plt.ylim(-0.1, 0.1)
    plt.legend(fontsize=15)
    plt.tight_layout()
    plt.savefig("eigen.png", dpi=200)
    plt.close()
