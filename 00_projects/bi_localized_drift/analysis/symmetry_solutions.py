from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    eq = 'PDNLS'
    disc = "aD:"
    initial_dir_data = str(disc) + 'mnustes_science/simulation_data/FD'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    x_grid = np.loadtxt(directory + '/X.txt', delimiter=',')
    #t_grid = np.loadtxt(directory + '/T.txt', delimiter=',')
    ZR = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    ZI = np.loadtxt(directory+ '/field_img.txt', delimiter=',')
    params = np.loadtxt(directory+ '/parameters.txt', delimiter=',')

    Z1 = ZR[-1, :] + 1j * ZI[-1, :]
    Z2 = np.flip(Z1) * np.exp(1j * np.pi / 2)
    Z3 = np.flip(Z2) * np.exp(1j * np.pi / 2)
    Z4 = np.flip(Z3) * np.exp(1j * np.pi / 2)

    t_rate = 100
    T = 3100
    dt = 0.012
    [tmin, tmax, dt] = [0, T, dt]
    t_grid = np.arange(tmin, tmax + dt, dt)
    dx = x_grid[1] - x_grid[0]
    xmin, xmax = x_grid[0], x_grid[-1]

    Nt = t_grid.shape[0]
    Nx = x_grid.shape[0]
    ies = [Z1, Z2, Z3, Z4]
    [alpha, beta, gamma_0, dist, sigma, mu, nu] = params
    print("N° of simulations: " + str(len(ies)))
    n = 0
    for i in ies:
        n = n + 1
        U_1_init = np.real(i)  #
        U_2_init = np.imag(i)  #

        [alpha_str, beta_str, mu_str, nu_str, sigma_str, gamma_str] = pdnlS_str_parameters([alpha, beta, mu, nu, sigma, gamma_0], 0)
        gamma_str = str(int(gamma_0 * 1000) * 0.001)
        nu_str = str(int(nu * 1000) * 0.001)
        mu_str = str(int(mu * 1000) * 0.001)
        print('gamma = ' + gamma_str)
        print('nu = ' + nu_str)
        print('mu = ' + mu_str)

        # Empaquetamiento de parametros, campos y derivadas para integración
        L = xmax - xmin
        D2 = sparse_DD_neumann(Nx, dx)
        operators = np.array([D2])
        fields_init = [U_1_init, U_2_init]
        grids = [t_grid, x_grid, 0]
        phi = np.pi
        gamma_complex = gamma_0 * (np.exp(- (x_grid - dist / 2) ** 2 / (2 * sigma ** 2)) + np.exp(1j * phi) * np.exp(- (x_grid + dist / 2) ** 2 / (2 * sigma ** 2)))
        gamma_real = np.real(gamma_complex)
        gamma_img = np.imag(gamma_complex)
        gamma = [gamma_real, gamma_img]
        parameters = [alpha, beta, gamma, mu, nu]

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
        U_complex = U1_light + 1j * U2_light
        t_light = time_grid

        # Definiendo variables finales
        modulo_light_1 = np.absolute(U_complex)
        arg_light_1 = np.angle(U_complex)
        arg_light_1 = (2 * np.pi + arg_light_1) * (arg_light_1 < 0) + arg_light_1 * (arg_light_1 > 0)
        analytical_signal_1 = hilbert(U1_light[-1, :])
        amplitude_envelope_1 = np.abs(analytical_signal_1)

        # Guardando datos
        subfile = pdnlS_bigauss_name([alpha, beta, mu, nu, sigma, gamma_0, dist], "ABMNSGD")
        parameters_np = np.array([alpha, beta, gamma_0, dist, sigma,  mu, nu])

        file = "/symmetry/" + str(n)
        if not os.path.exists(directory + file):
            os.makedirs(directory + file)
        np.savetxt(directory + file + '/field_real.txt', U1_light, delimiter=',')
        np.savetxt(directory + file + '/field_img.txt', U2_light, delimiter=',')
        np.savetxt(directory + file + '/parameters.txt', parameters_np, delimiter=',')
        np.savetxt(directory + file + '/X.txt', x_grid, delimiter=',')
        np.savetxt(directory + file + '/T.txt', t_light, delimiter=',')

        # Guardando gráficos
        pcm = plt.pcolormesh(x_grid, t_light, modulo_light_1, cmap=parula_map, shading='auto')
        cbar = plt.colorbar(pcm, shrink=1)
        cbar.set_label('$|A|$', rotation=0, size=20, labelpad=-27, y=1.1)
        plt.xlim([x_grid[0], x_grid[-1]])
        plt.xlabel('$x$', size='20')
        plt.ylabel('$t$', size='20')
        plt.grid(linestyle='--', alpha=0.5)
        plt.savefig(directory + file + '/module_spacetime.png', dpi=300)
        plt.close()

        pcm = plt.pcolormesh(x_grid, t_light, arg_light_1, cmap=parula_map, shading='auto')
        cbar = plt.colorbar(pcm, shrink=1)
        cbar.set_label('$\\textrm{arg}(A)$', rotation=0, size=20, labelpad=-20, y=1.1)
        plt.xlim([x_grid[0], x_grid[-1]])
        plt.xlabel('$x$', size='20')
        plt.ylabel('$t$', size='20')
        plt.grid(linestyle='--', alpha=0.5)
        plt.savefig(directory + file + '/arg_spacetime.png', dpi=300)
        plt.close()