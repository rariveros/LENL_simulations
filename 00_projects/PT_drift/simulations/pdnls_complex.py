import matplotlib.pyplot as plt

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


    [alpha, beta, gamma_0, mu, nu, sigma, phi] = params_00
    gamma = -(gaussian(x_grid, gamma_0, -d / 2, sigma) - gaussian(x_grid, gamma_0, d / 2, sigma))

    plt.plot(x_grid, np.real(PHI_R), color="b", label="$\phi_{R}^{R}$")
    plt.plot(x_grid, np.imag(PHI_R), color="r", label="$\phi_{R}^{I}$")
    #plt.plot(x_grid, np.real(PHI_L), color="b", linestyle="--", label="$\phi_{L}^{R}$")
    #plt.plot(x_grid, np.imag(PHI_L), color="r", linestyle="--", label="$\phi_{L}^{I}$")
    plt.plot(x_grid, gamma, color="k", label="$\gamma(x)$")

    plt.legend(fontsize=15)
    plt.show()
    plt.close()

    distances = np.arange(50, 150, 2)
    PHIs = []
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
        PHI_L = 1j * PHI_L
        PHI_R = 1j * PHI_R
        PHI = PHI_L + PHI_R
        phi_01 = PHI_L - 1j * PHI_R
        phi_02 = PHI_L + 1j * PHI_R

        [alpha, beta, gamma_0, mu, nu, sigma, phi] = params_00
        gamma = -(gaussian(x_grid, gamma_0, -d / 2, sigma) - gaussian(x_grid, gamma_0, d / 2, sigma))

        Gamma_0 = integrate.simpson(phi_01 * phi_02, x_grid)
        Gamma_1 = integrate.simpson(phi_01 ** 2, x_grid)
        Gamma_2 = integrate.simpson(phi_02 ** 2, x_grid)

        ddx_phi_01 = np.append(0, np.append(np.diff(np.diff(phi_01)), 0))
        ddx_phi_02 = np.append(0, np.append(np.diff(np.diff(phi_02)), 0))

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

        sigma = [Sigma_11, Sigma_12, Sigma_21, Sigma_22]
        pi = [Pi_11, Pi_12, Pi_21, Pi_22]
        delta = [Delta_11, Delta_21, Delta_31, Delta_41, Delta_51, Delta_61, Delta_12, Delta_22, Delta_32, Delta_42, Delta_52, Delta_62]
        parameters = [alpha, beta, mu, nu, sigma, pi, delta]

        PHIs.append(PHI)

        # Definición de la grilla
        [tmin, tmax, dt] = [0, T, dt]
        t_grid = np.arange(tmin, tmax + dt, dt)
        T = tmax
        Nt = t_grid.shape[0]

        # Initial Conditions Pattern
        U_1_init = 0.01 * (np.random.rand(1) + 1j * np.random.rand(1))
        U_2_init = 0.01 * (np.random.rand(1) + 1j * np.random.rand(1))
        operators = [0]

        # Empaquetamiento de parametros, campos y derivadas para integración
        fields_init = [U_1_init, U_2_init]
        grids = [t_grid, 0, 0]

        # Midiendo tiempo inicial
        now = datetime.datetime.now()
        print('Hora de Inicio: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
        time_init = time.time()

        final_fields, fields_history, time_grid = RK4_FD(eq, fields_init, parameters, grids, dt, Nt, operators, t_rate) #INTEGRACION EN EL TIEMPO

        now = datetime.datetime.now()
        print('Hora de Término: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
        time_fin = time.time()
        print(str(time_fin - time_init) + ' seg')

        # Reobteniendo campos
        U1_light = np.array(fields_history)[:, 0]
        U2_light = np.array(fields_history)[:, 1]
        t_light = time_grid

        UL = (U1_light + U2_light) * ((1 - 1j) / 2)
        UR = - 1j * (U1_light - U2_light) * ((1 + 1j) / 2)

        if d == 28:
            Z_recons = []
            for i in range(len(t_light)):
                Z_recons_i = PHI_R * UR[i] + PHI_L * UL[i]
                Z_recons.append(Z_recons_i)
            pcm = plt.pcolormesh(x_grid, t_light, np.abs(Z_recons), cmap=parula_map, shading='auto')
            cbar = plt.colorbar(pcm, shrink=1)
            cbar.set_label('$|A|$', rotation=0, size=20, labelpad=-27, y=1.1)
            plt.xlim([x_grid[0], x_grid[-1]])
            plt.xlabel('$x$', size='20')
            plt.ylabel('$t$', size='20')
            plt.grid(linestyle='--', alpha=0.5)
            plt.show()
            plt.close()


        save_directory = directory + "/analysis"

        if not os.path.exists(save_directory):
            os.makedirs(save_directory)

        plt.plot(t_light, np.real(UL), c="b", label="$\\textrm{Re}(\psi_L)$")
        plt.plot(t_light, np.imag(UL), c="r", label="$\\textrm{Im}(\psi_L)$")
        plt.plot(t_light, np.real(UR), c="b", linestyle="--", label="$\\textrm{Re}(\psi_R)$")
        plt.plot(t_light, np.imag(UR), c="r", linestyle="--", label="$\\textrm{Im}(\psi_R)$")
        #plt.plot(t_light, np.abs(UL), c="green", label="$|\psi_L|$")
        #plt.plot(t_light, np.abs(UR), c="purple", label="$|\psi_R|$")
        plt.xlabel("$t$")
        plt.legend()
        plt.savefig(save_directory + "/d=" + str(d) + ".png")
        #plt.show()
        plt.close()

    fig, ((ax1a, ax2a, ax3a), (ax1b, ax2b, ax3b)) = plt.subplots(nrows=2, ncols=3)

    pc_01b = ax1b.pcolor(x_grid, distances, np.real(PHIs), cmap=parula_map)
    ax1b.set_ylabel("$d$", fontsize=15)
    ax1b.tick_params(axis="x", direction="in", labeltop=True, labelbottom=False, top=True, bottom=True)
    ax1b.tick_params(axis="y", direction="in", left=True, right=True)
    ax1b.grid(alpha=0.2, color="k")
    ax1b.set_xlim(-150, 150)

    pc_02b = ax2b.pcolor(x_grid, distances, np.imag(PHIs), cmap=parula_map)
    ax2b.set_yticklabels([])
    ax2b.tick_params(axis="x", direction="in", labeltop=True, labelbottom=False, top=True, bottom=True)
    ax2b.tick_params(axis="y", direction="in", left=True, right=True)
    ax2b.grid(alpha=0.2, color="k")
    ax2b.set_xlim(-150, 150)

    pc_03b = ax3b.pcolor(x_grid, distances, np.abs(PHIs), cmap=parula_map)
    ax3b.set_yticklabels([])
    ax3b.tick_params(axis="x", direction="in", labeltop=True, labelbottom=False, top=True, bottom=True)
    ax3b.tick_params(axis="y", direction="in", left=True, right=True)
    ax3b.grid(alpha=0.2, color="k")
    ax3b.set_xlim(-150, 150)

    plt.show()
    plt.close()