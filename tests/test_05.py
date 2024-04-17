from functions import *
from back_process import *
from time_integrators import *


if __name__ == '__main__':
    disco = 'D:/'
    initial_dir_data = str(disco) + 'Users/mnustes_science/PT_fluids/mnustes_science/simulation_data'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    x_grid = np.loadtxt(directory + '/X.txt', delimiter=',')
    params_00 = np.loadtxt(directory + '/parameters.txt', delimiter=',')

    [alpha, beta, gamma_0, mu, nu, sigma_i] = params_00
    dx = x_grid[1] - x_grid[0]
    gamma = gamma_0 #* np.exp(- x_grid ** 2 / (2 * sigma_i ** 2))
    mu = 0.1
    A = np.arange(-0.9, 5, 0.01)
    F = []
    for i in range(len(A)):
        phi_01_i = Z_r[-1, :] / (1 + A[i])
        phi_02_i = Z_i[-1, :] / (1 + A[i])

        Dphi_01_i = np.append(np.diff(phi_01_i) / (dx), 0)
        Dphi_02_i = np.append(np.diff(phi_02_i) / (dx), 0)

        f_ij = (-(alpha / 2) * (Dphi_01_i ** 2 + Dphi_02_i ** 2) + (beta / 4) * (phi_01_i ** 2 + phi_02_i ** 2) ** 2 + (nu / 2) * (phi_01_i ** 2 + phi_02_i ** 2) + gamma * phi_01_i * phi_02_i) / (2 * mu)
        #if A[i] == 1:
        #    plt.plot(x_grid, f_ij)
        #    plt.show()
        #    plt.close()
        F_i = integrate.simpson(f_ij, x_grid)
        F.append(F_i)
    F = np.array(F)

    plt.plot(A, np.append(np.diff(F), 0))
    plt.show()