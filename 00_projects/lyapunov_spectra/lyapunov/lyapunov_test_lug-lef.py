from functions import *
from back_process import *
from time_integrators import *
from jacobians import *
from directories_lyap import *

if __name__ == '__main__':
    eq = "lug_lef"
    disco = 'E:/'
    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FD' + main_directory
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    directories = os.listdir(working_directory)

    # INITIAL TIME
    now = datetime.datetime.now()
    print('Hora de Inicio: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
    time_init = time.time()

    for directory_i in directories:
        directory = working_directory + "/" + directory_i
        Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
        Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
        x_grid = np.loadtxt(directory + '/X.txt', delimiter=',')
        t_grid = np.loadtxt(directory + '/T.txt', delimiter=',')

        # DEFINING TEMPORAL DOMAIN TO ANALYZE
        Nx = len(x_grid)
        Nt = len(t_grid)
        dx = x_grid[1] - x_grid[0]
        dt = 0.005
        ti = int(0.6 * Nt)
        tf = -1
        t_grid = t_grid[ti:tf]
        Nt = len(t_grid)
        Z_r = Z_r[ti:tf, :]
        Z_i = Z_i[ti:tf, :]

        # DEFINING INITIAL CONDITIONS
        print("dt = " + str(dt))
        lyap = []
        N_condit = int(2 * Nx)
        print("N_c = " + str(N_condit))
        U_init = np.random.rand(2 * Nx, N_condit) - 0.5
        Q, R = np.linalg.qr(U_init)

        # DEFINING PARAMETERS
        parameters = np.loadtxt(directory + '/parameters.txt', delimiter=',')
        DD = sparse_DD(Nx, dx)
        subfile = "/det=" + str(parameters[0]) + "/E=" + str(parameters[1])
        save_directory = initial_dir_data + subdir_exponent + "/test_01" + subfile

        I = np.eye(2 * Nx)

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