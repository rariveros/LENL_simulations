from functions import *
from back_process import *
from time_integrators import *
from directories_lyap import *

if __name__ == '__main__':
    disco = 'E:/'
    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FD' + main_directory + directoy_sim + subdir_branching
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    directories = [name for name in os.listdir(working_directory) if os.path.isdir(os.path.join(working_directory, name))]
    freq = []
    ampd = []
    av_power = []
    param = []
    for directory_i in directories:
        directory = working_directory + "/" + directory_i
        print("-----------" + directory_i + "-----------")
        Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
        Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
        X = np.loadtxt(directory + '/X.txt', delimiter=',')
        T = np.loadtxt(directory + '/T.txt', delimiter=',')

        Nx = len(X)
        Nt = len(T)
        dt = T[1] - T[0]
        dx = X[1] - X[0]
        L = X[-1] - X[0]
        sigma_i = 60
        mu = 0.1
        nu = 0.056
        alpha = 5.721
        sigma = sigma_i #2.3548 * sigma_i#np.sqrt(np.sqrt(nu * alpha) * (sigma_i / mu))
        sigma_window_l = int((Nx - (sigma / L) * Nx) / 2)
        sigma_window_r = int(sigma_window_l + (sigma / L) * Nx)

        # DEFINE REGION OF INTEREST +- SIGMA/2

        ti = int(0.5 * Nt)
        tf = -1
        T = T[ti:tf]
        Nt = len(T)
        Z_r = Z_r[ti:tf, :]
        Z_i = Z_i[ti:tf, :]

        Z_complex = Z_r + 1j * Z_i

        # Definiendo variables finales
        Z_modulo= np.absolute(Z_complex)
        Z_arg = np.angle(Z_complex)
        Z_arg = (2 * np.pi + Z_arg) * (Z_arg < 0) + Z_arg * (Z_arg > 0)

        D = sparse_D_neumann(Nx, dx)
        DD = sparse_DD_neumann(Nx, dx)

        D1_Z = Der(D, np.transpose(Z_modulo))
        D2_Z = Der(DD, np.transpose(Z_modulo))
        D1_Z = np.transpose(D1_Z)
        D2_Z = np.transpose(D2_Z)

        J_points_init = []
        K = 0
        for j in range(sigma_window_l, sigma_window_r):
            if np.sign(D1_Z[0, j]) != np.sign(D1_Z[0, j - 1]) and D2_Z[0, j] < 0 and j != 0:
                K = K + 1
                J_points_init.append(j)
        print(K)

        Z_r_points = []
        Z_i_points = []
        amplitude_points = []
        x_points = []
        t_points = []
        window_size = 4
        for k in range(K):
            window_left = J_points_init[k] - int(window_size / 2)
            window_right = J_points_init[k] + int(window_size / 2)
            x_points_k = []
            amplitude_points_k = []
            t_points_k = []
            Z_r_points_k = []
            Z_i_points_k = []
            for i in range(Nt):
                for j in range(window_size):
                    if np.sign(D1_Z[i, window_left + j]) != np.sign(D1_Z[i, window_left + j - 1]) and D2_Z[i, window_left + j] < 0 and j != 0:
                        Z_r_points_k.append(Z_r[i, window_left + j])
                        Z_i_points_k.append(Z_i[i, window_left + j])
                        amplitude_points_k.append(Z_modulo[i, window_left + j])
                        x_points_k.append(X[window_left + j])
                        t_points_k.append(T[i])
                        window_left = window_left + j - int(window_size / 2)
                        window_right = window_left + j + int(window_size / 2)
            Z_r_points.append(Z_r_points_k)
            Z_i_points.append(Z_i_points_k)
            amplitude_points.append(amplitude_points_k)
            t_points.append(t_points_k)
            x_points.append(x_points_k)
        n = int(K / 2)
        Z_r_points_np = np.array(Z_r_points[n][:])
        Z_i_points_np = np.array(Z_i_points[n][:])
        x_points_np = np.array(x_points[n][:])
        t_points_np = np.array(t_points[n][:])
        np.savetxt(directory + '/Z_r_points.txt', Z_r_points_np, delimiter=',')
        np.savetxt(directory + '/Z_i_points.txt', Z_i_points_np, delimiter=',')
        np.savetxt(directory + '/x_points.txt', x_points_np, delimiter=',')
        np.savetxt(directory + '/t_points.txt', t_points_np, delimiter=',')


        Nt = len(t_points_np)
        dt = np.abs(t_points_np[1] - t_points_np[0])

        av_power_i = np.sum(np.abs((Z_r_points_np - np.mean(Z_r_points_np)) * (Z_r_points_np - np.mean(Z_r_points_np))) ** 2) / Nt

        yf_i = fft(Z_r_points_np)
        yf_i = 2.0 / Nt * np.abs(yf_i[0:Nt // 2])
        xf_i = fftfreq(Nt, dt)[:Nt // 2]

        yf_i = yf_i[1:]
        xf_i = xf_i[1:]

        freq_i = xf_i[np.argmax(yf_i)]
        ampd_i = yf_i[np.argmax(yf_i)]

        plt.plot(xf_i, np.log(yf_i), color="k")
        plt.xlim(0, 1)
        plt.xlabel("$\Omega$", fontsize=20)
        plt.ylabel("$\\textrm{Power Spectrum}$", fontsize=20)
        plt.grid(alpha=0.5)
        plt.tight_layout()
        plt.savefig(directory + '/fft_amplitude.png', dpi=300)
        plt.close()


        freq.append(freq_i)
        ampd.append(ampd_i)
        av_power.append(av_power_i)
        print(freq)
        param.append(float(directory_i.split("=")[-1]))
    np.savetxt(working_directory + '/param.txt', param, delimiter=',')
    np.savetxt(working_directory + '/freq.txt', freq, delimiter=',')
    np.savetxt(working_directory + '/ampd.txt', ampd, delimiter=',')
    np.savetxt(working_directory + '/power.txt', av_power, delimiter=',')

    plt.scatter(param, freq, c="k", zorder=10)
    plt.grid(alpha=0.5, zorder=1)
    plt.ylim([0, 1.1 * np.amax(freq)])
    plt.xlabel("$\gamma_0$", fontsize=20)
    plt.ylabel("$f$", fontsize=20)
    plt.tight_layout()
    plt.savefig(working_directory + '/freq_oscillations.png', dpi=300)
    plt.close()

    plt.scatter(param, ampd, c="k", zorder=10)
    plt.grid(alpha=0.5, zorder=1)
    plt.ylim([0, 1.1 * np.amax(ampd)])
    plt.xlabel("$\gamma_0$", fontsize=20)
    plt.ylabel("$\delta A_R$", fontsize=20)
    plt.tight_layout()
    plt.savefig(working_directory + '/ampd_oscillations.png', dpi=300)
    plt.close()

    plt.scatter(param, av_power, c="k", zorder=10)
    plt.grid(alpha=0.5, zorder=1)
    plt.ylim([0, 1.1 * np.amax(av_power)])
    plt.xlabel("$\gamma_0$", fontsize=20)
    plt.ylabel("$\\textrm{Average Power}$", fontsize=20)
    plt.tight_layout()
    plt.savefig(working_directory + '/power_oscillations.png', dpi=300)
    plt.close()

