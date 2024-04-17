from functions import *
from back_process import *
from time_integrators import *

def gaussian(x, A, x_0, sigma):
    return A * np.exp(- (x - x_0) ** 2 / (2 * sigma ** 2))

if __name__ == '__main__':
    disco = 'D:/'
    initial_dir_data = str(disco) + 'Users/mnustes_science/PT_fluids/mnustes_science/simulation_data'
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    directories = [name for name in os.listdir(working_directory) if os.path.isdir(os.path.join(working_directory, name))]
    Z_Rs = []
    Z_Is = []
    Z_mods = []
    dists = []

    ##### ANSATZ #####
    Z_r_00 = np.loadtxt('E:/mnustes_science/simulation_data/FD/PT_drift/responses/alpha=6.524/beta=1.000/mu=0.100/nu=0.018/sigma=6.000/gamma=0.185/phi00/field_real_0.txt', delimiter=',')
    Z_i_00 = np.loadtxt('E:/mnustes_science/simulation_data/FD/PT_drift/responses/alpha=6.524/beta=1.000/mu=0.100/nu=0.018/sigma=6.000/gamma=0.185/phi00/field_img_0.txt', delimiter=',')
    params_00 = np.loadtxt('E:/mnustes_science/simulation_data/FD/PT_drift/responses/alpha=6.524/beta=1.000/mu=0.100/nu=0.018/sigma=6.000/gamma=0.185/phi00/parameters.txt', delimiter=',')

    Z_r_01 = np.loadtxt('E:/mnustes_science/simulation_data/FD/PT_drift/responses/alpha=6.524/beta=1.000/mu=0.100/nu=0.018/sigma=6.000/gamma=0.185/phi01/field_real_0.txt', delimiter=',')
    Z_i_01 = np.loadtxt('E:/mnustes_science/simulation_data/FD/PT_drift/responses/alpha=6.524/beta=1.000/mu=0.100/nu=0.018/sigma=6.000/gamma=0.185/phi01/field_img_0.txt', delimiter=',')
    x_grid_phi = np.loadtxt('E:/mnustes_science/simulation_data/FD/PT_drift/responses/alpha=6.524/beta=1.000/mu=0.100/nu=0.018/sigma=6.000/gamma=0.185/phi01/X.txt', delimiter=',')
    params_01 = np.loadtxt('E:/mnustes_science/simulation_data/FD/PT_drift/responses/alpha=6.524/beta=1.000/mu=0.100/nu=0.018/sigma=6.000/gamma=0.185/phi01/parameters.txt', delimiter=',')

    PHI_Rs = []
    PHI_Is = []
    PHI_mods = []
    for directory in directories:
        print(directory)
        working_directory_i = working_directory + "/" + directory
        Z_r = np.loadtxt(working_directory_i + '/field_real.txt', delimiter=',')
        Z_i = np.loadtxt(working_directory_i + '/field_img.txt', delimiter=',')
        x_grid = np.loadtxt(working_directory_i + '/X.txt', delimiter=',')
        params = np.loadtxt(working_directory_i + '/parameters.txt', delimiter=',')

        d = params[3]
        Z_r_final = Z_r[-1, :]
        Z_i_final = Z_i[-1, :]
        Z_mod_final = np.abs(Z_r_final + 1j * Z_i_final)

        X_L = - d / 2
        X_R = + d / 2
        J_L = np.argmin(np.abs(x_grid_phi - X_L))
        J_R = np.argmin(np.abs(x_grid_phi - X_R))
        J_center = np.argmin(np.abs(x_grid_phi))
        Delta_J_L = J_center - J_L
        Delta_J_R = J_R - J_center
        PHI_L = np.append(Z_r_00[Delta_J_L:], np.zeros(Delta_J_L)) + 1j * np.append(Z_i_00[Delta_J_L:], np.zeros(Delta_J_L))
        PHI_R = np.append(np.zeros(Delta_J_R), Z_r_01[:-Delta_J_R]) + 1j * np.append(np.zeros(Delta_J_R), Z_i_01[:-Delta_J_R])
        PHI = -1j * (PHI_L + PHI_R)

        Z_Rs.append(Z_r_final)
        Z_Is.append(Z_i_final)
        Z_mods.append(Z_mod_final)
        dists.append(d)

        PHI_Rs.append(np.real(PHI))
        PHI_Is.append(np.imag(PHI))
        PHI_mods.append(np.abs(PHI))

        if d == 149:
            fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
            ax1.plot(x_grid, Z_r_final, label="$\psi_R$")
            ax1.plot(x_grid_phi, np.real(PHI), label="$\phi_R$")
            ax1.set_xlim([-150, 150])
            ax1.legend()

            ax2.plot(x_grid, Z_i_final, label="$\psi_I$")
            ax2.plot(x_grid_phi, np.imag(PHI), label="$\phi_I$")
            ax2.set_xlim([-150, 150])
            ax2.legend()

            plt.show()
            plt.close()

    sorted_indices = np.argsort(dists)
    Z_Rs_sorted = []
    Z_Is_sorted = []
    Z_mods_sorted = []
    PHI_Rs_sorted = []
    PHI_Is_sorted = []
    PHI_mods_sorted = []
    dists_sorted = []

    for i in range(0, len(sorted_indices)):
        Z_Rs_sorted.append(Z_Rs[sorted_indices[i]])
        Z_Is_sorted.append(Z_Is[sorted_indices[i]])
        Z_mods_sorted.append(Z_mods[sorted_indices[i]])
        PHI_Rs_sorted.append(PHI_Rs[sorted_indices[i]])
        PHI_Is_sorted.append(PHI_Is[sorted_indices[i]])
        PHI_mods_sorted.append(PHI_mods[sorted_indices[i]])
        dists_sorted.append(dists[sorted_indices[i]])
    x_L = x_grid[0]
    x_R = x_grid[-1]
    fig, ((ax1a, ax2a, ax3a), (ax1b, ax2b, ax3b)) = plt.subplots(nrows=2, ncols=3)

    pc_01a = ax1a.pcolor(x_grid, dists, Z_Rs_sorted, cmap=parula_map)
    ax1a.set_ylabel("$d$", fontsize=15)
    ax1a.tick_params(axis="x", direction="in", labeltop=True, labelbottom=False, top=True, bottom=True)
    ax1a.tick_params(axis="y", direction="in", left=True, right=True)
    ax1a.grid(alpha=0.2, color="k")

    pc_02a = ax2a.pcolor(x_grid, dists, Z_Is_sorted, cmap=parula_map)
    ax2a.set_yticklabels([])
    ax2a.tick_params(axis="x", direction="in", labeltop=True, labelbottom=False, top=True, bottom=True)
    ax2a.tick_params(axis="y", direction="in", left=True, right=True)
    ax2a.grid(alpha=0.2, color="k")

    pc_03a = ax3a.pcolor(x_grid, dists, Z_mods_sorted, cmap=parula_map)
    ax3a.set_yticklabels([])
    ax3a.tick_params(axis="x", direction="in", labeltop=True, labelbottom=False, top=True, bottom=True)
    ax3a.tick_params(axis="y", direction="in", left=True, right=True)
    ax3a.grid(alpha=0.2, color="k")

    #####

    pc_01b = ax1b.pcolor(x_grid_phi, dists, PHI_Rs_sorted, cmap=parula_map)
    ax1b.set_ylabel("$d$", fontsize=15)
    ax1b.tick_params(axis="x", direction="in", labeltop=True, labelbottom=False, top=True, bottom=True)
    ax1b.tick_params(axis="y", direction="in", left=True, right=True)
    ax1b.grid(alpha=0.2, color="k")
    ax1b.set_xlim(x_L, x_R)

    b = ax2b.pcolor(x_grid_phi, dists, PHI_Is_sorted, cmap=parula_map)
    ax2b.set_yticklabels([])
    ax2b.tick_params(axis="x", direction="in", labeltop=True, labelbottom=False, top=True, bottom=True)
    ax2b.tick_params(axis="y", direction="in", left=True, right=True)
    ax2b.grid(alpha=0.2, color="k")
    ax2b.set_xlim(x_L, x_R)

    pc_03b = ax3b.pcolor(x_grid_phi, dists, PHI_mods_sorted, cmap=parula_map)
    ax3b.set_yticklabels([])
    ax3b.tick_params(axis="x", direction="in", labeltop=True, labelbottom=False, top=True, bottom=True)
    ax3b.tick_params(axis="y", direction="in", left=True, right=True)
    ax3b.grid(alpha=0.2, color="k")
    ax3b.set_xlim(x_L, x_R)

    plt.show()
    plt.close()