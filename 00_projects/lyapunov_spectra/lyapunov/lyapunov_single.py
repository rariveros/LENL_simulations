from functions import *
from back_process import *
from time_integrators import *
from jacobians import *
from directories_lyap import *

if __name__ == '__main__':
    eq = 'pdnlS'
    disco = 'E:/'
    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FD' + main_directory
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    x_grid = np.loadtxt(directory + '/X.txt', delimiter=',')
    t_grid = np.loadtxt(directory + '/T.txt', delimiter=',')

    # DEFINING TEMPORAL DOMAIN TO ANALYZE
    Nx = len(x_grid)
    Nt = len(t_grid)
    dx = x_grid[1] - x_grid[0]
    dt = 0.02
    ti = int(0.8 * Nt)
    tf = -1
    t_grid = t_grid[ti:tf]
    Nt = len(t_grid)
    Z_r = Z_r[ti:tf, :]
    Z_i = Z_i[ti:tf, :]

    for i in range(1):
        # DEFINING INITIAL CONDITIONS
        lyap = []
        N_condit = int(0.3 * 2 * Nx)
        print("dt = " + str(dt))
        print("N_c = " + str(N_condit))
        print("Nt = " + str(Nt))
        U_init = np.random.rand(2 * Nx, N_condit) - 0.5
        Q, R = np.linalg.qr(U_init)

        # DEFINING PARAMETERS
        sigma = float(directory.split("=")[-1])
        [alpha, beta, gamma_0, mu_0, nu] = np.loadtxt(directory + '/parameters.txt', delimiter=',')
        print("sigma = " + str(sigma))
        gamma_1 = gamma_0 * np.exp(- x_grid ** 2 / (2 * sigma ** 2))
        gamma_2 = 0
        gamma = [gamma_1, gamma_2]
        mu = mu_0 * np.ones(Nx)
        mu[0:10] = 10
        mu[-10:-1] = 10
        subfile = pdnlS_name([alpha, beta, mu_0, nu, sigma, gamma_0], "ABMNGS")
        parameters = [alpha, beta, gamma, mu, nu]
        save_directory = initial_dir_data + subdir_exponent + "/test_03" + subfile + "/0" + str(i)
        DD = sparse_DD_neumann(Nx, dx)
        I = np.eye(2 * Nx)

        # INITIAL TIME
        now = datetime.datetime.now()
        print('Hora de Inicio: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
        time_init = time.time()

        for i in range(Nt):
            U_1 = Z_r[i, :]
            U_2 = Z_i[i, :]
            J = jacobians_FD(eq, [U_1, U_2], t_grid, x_grid, [0], parameters, [DD])
            Q_new = time_propagator("I_Jdt_RK4", I, J, Q, dt)
            Q, R = np.linalg.qr(Q_new)
            Q = Q
            lyap.append(np.log(np.absolute(R.diagonal())) / dt)
        lyap = np.array(lyap)
        lyap_mean = []
        average_Nt_i = int(0.3 * Nt)
        for i in range(N_condit):
            lyap_mean.append(np.sum(lyap[average_Nt_i:-1, i]))
            plt.plot(lyap[:, i])
        lyap_mean = np.array(lyap_mean) / (Nt - average_Nt_i)
        np.savetxt(directory + "/lyap.txt", lyap_mean)
        plt.tight_layout()
        plt.savefig(directory + '/lyap_espectrum.png', dpi=300)
        plt.close()

        x = np.arange(N_condit)
        plt.plot(lyap_mean, label="$dt = " + str(np.round(dt, 3)) + "$", linewidth=0.5)
        plt.legend()
        plt.title("$\\textrm{Lyapunov Spectrum}$", size='15')
        plt.legend(fontsize=8)
        plt.xlabel('$\\textrm{N° Initial Conditions}$', size='15')
        plt.xticks(fontsize=15)
        plt.xlim([0, 1200])
        plt.ylabel('$\lambda$', size='20')
        plt.yticks(fontsize=15)
        plt.grid(alpha=0.2)
        plt.tight_layout()
        plt.savefig(directory + '/lyap_espectrums.png', dpi=300)
        plt.close()
        if not os.path.exists(save_directory):
            os.makedirs(save_directory)
        np.savetxt(save_directory + '/x_lyap_' + str(N_condit) + '.txt', x, delimiter=',')
        np.savetxt(save_directory + '/lyap_mean_' + str(N_condit) + '.txt', lyap_mean, delimiter=',')

        # FINAL TIME
        now = datetime.datetime.now()
        print('Hora de Término: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
        time_fin = time.time()
        print(str(np.around((time_fin - time_init) / 60, 3)) + ' min')
        print("-------------------------------------")