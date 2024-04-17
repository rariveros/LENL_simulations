from back_process import *
from jacobians import *

def gaussian(x, A, x_0, sigma):
    return A * np.exp(- (x - x_0) ** 2 / (2 * sigma ** 2))
def resample(old_arrays, old_x, N_resample):
    new_arrays = []
    for i in range(len(old_arrays)):
        [new_array_i, new_x] = signal.resample(old_arrays[i], int(len(old_x) * (N_resample)), old_x)
        new_arrays.append(new_array_i)
    return new_arrays, new_x
if __name__ == '__main__':
    disco = 'D:/'
    eq = 'pdnlS'

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
    #N_resample = 10
    #[Z_r_00, Z_i_00, Z_r_01, Z_i_01], x_grid = resample([Z_r_00, Z_i_00, Z_r_01, Z_i_01], x_grid, N_resample)
    Nx = len(x_grid)
    [alpha, beta, gamma_0, mu, nu, sigma, phi] = params_00
    distances = np.arange(30.4, 31.1, 0.1)
    colors = np.arange(0, 1 + 1 / (len(distances) - 1), 1 / (len(distances) + 1))
    for i in range(len(distances)):
        d = distances[i]
        print(d)
        gamma_real = -(gaussian(x_grid, gamma_0, -d / 2, sigma) - gaussian(x_grid, gamma_0, d / 2, sigma))
        gamma_img = 0
        gamma = [gamma_real, gamma_img]
        parameters = [alpha, beta, gamma, mu, nu]
        DD = sparse_DD_neumann(Nx, dx)
        X_L = - d / 2
        X_R = + d / 2
        J_L = np.argmin(np.abs(x_grid - X_L))
        J_R = np.argmin(np.abs(x_grid - X_R))
        J_center = np.argmin(np.abs(x_grid))
        Delta_J_L = J_center - J_L
        Delta_J_R = J_R - J_center
        PHI_L = np.append(Z_r_00[Delta_J_L:], np.zeros(Delta_J_L)) + 1j * np.append(Z_i_00[Delta_J_L:], np.zeros(Delta_J_L))
        PHI_R = np.append(np.zeros(Delta_J_R), Z_r_01[:-Delta_J_R]) + 1j * np.append(np.zeros(Delta_J_R),
                                                                                     Z_i_01[:-Delta_J_R])
        #PHI_L = 1j * PHI_L
        #PHI_R = 1j * PHI_R

        PHI = np.flip((PHI_L + PHI_R))
        #plt.plot(x_grid, np.real(PHI))
        #plt.plot(x_grid, np.imag(PHI))
        #plt.show()
        #plt.close()

        J = jacobians_FD(eq, [np.real(PHI), np.imag(PHI)], [0], x_grid, [0], parameters, [DD])
        eigenvalues, eigenvectors = np.linalg.eig(J)
        arg_max = np.argmax(np.real(eigenvalues))
        plt.scatter(np.real(eigenvalues[arg_max]), np.imag(eigenvalues[arg_max]), marker="D", c=(colors[i], 0, 1 - colors[i]), zorder=10)
    plt.grid(alpha=0.4, zorder=0)
    # plt.xlim([-0.1, 0.1])
    # plt.ylim([-1.2, 1.2])
    plt.xlabel("$\\textrm{Re}\{\lambda_i\}$", fontsize=18)
    plt.ylabel("$\\textrm{Im}\{\lambda_i\}$", fontsize=18)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.vlines(0, -1.2, 1.2, colors="k")
    plt.hlines(0, -0.1, 0.1, colors="k")
    plt.tight_layout()
    plt.show()
    plt.close()