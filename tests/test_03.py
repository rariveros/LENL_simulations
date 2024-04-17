from functions import *
from back_process import *
from time_integrators import *

def gaussian(x, A, x_0, sigma):
    return A * np.exp(- (x - x_0) ** 2 / (2 * sigma ** 2))

if __name__ == '__main__':

    eq = 'PT_dimer'
    t_rate = 1
    dt = 1
    T = 2000

    disco = 'D:/'
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

    d = 100
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



    plt.plot(x_grid, PHI_R)
    plt.plot(x_grid_resampled, PHI_R_resampled)
    plt.scatter(x_grid, PHI_R, s=20)
    plt.scatter(x_grid_resampled, PHI_R_resampled, s=10)
    plt.show()