from functions import *
from back_process import *
from time_integrators import *


if __name__ == '__main__':
    disco = 'E:/'
    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FD'
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    directories = [name for name in os.listdir(working_directory) if os.path.isdir(os.path.join(working_directory, name))]
    process_directory = r"D:/mnustes_science/simulation_data/FD/PDNLS_oscillatory/alpha=5.721/beta=1.000/mu=0.100/process/nu=0.027"
    for directory_i in directories:
        directory = working_directory + "/" + directory_i
        sub_directories = [name for name in os.listdir(directory) if os.path.isdir(os.path.join(directory, name))]
        print("-----------" + directory_i + "-----------")
        yf = []
        xf = []
        freq = []
        ampd = []
        av_power = []
        param = []
        sigma_i = float(directory_i.split("=")[-1])
        for directory_j in sub_directories:
            subdirectory = directory + "/" + directory_j
            Z_r_points_np = np.loadtxt(subdirectory + '/Z_r_points.txt', delimiter=',')
            x_points_np = np.loadtxt(subdirectory + '/x_points.txt', delimiter=',')
            t_points_np = np.loadtxt(subdirectory + '/t_points.txt', delimiter=',')
            Nt = len(t_points_np)
            dt = np.abs(t_points_np[1] - t_points_np[0])

            av_power_i = np.sum(
                np.abs((Z_r_points_np - np.mean(Z_r_points_np)) * (Z_r_points_np - np.mean(Z_r_points_np))) ** 2) / Nt

            yf_i = fft(Z_r_points_np)
            yf_i = 2.0 / Nt * np.abs(yf_i[0:Nt // 2])
            xf_i = fftfreq(Nt, dt)[:Nt // 2]

            yf_i_paramax = yf_i[50:]
            xf_i_paramax = xf_i[50:]

            freq_i = xf_i[np.argmax(yf_i_paramax)]
            ampd_i = yf_i[np.argmax(yf_i_paramax)]
            freq.append(freq_i)
            ampd.append(ampd_i)
            av_power.append(av_power_i)

            np.savetxt(subdirectory + '/yfft.txt', yf_i, delimiter=',')
            np.savetxt(subdirectory + '/xfft.txt', xf_i, delimiter=',')

            plt.plot(xf_i, np.log(yf_i), color="k")
            plt.xlim(0, 1)
            plt.xlabel("$\Omega$", fontsize=20)
            plt.ylabel("$\\textrm{Power Spectrum}$", fontsize=20)
            plt.grid(alpha=0.5)
            plt.tight_layout()
            plt.savefig(subdirectory + '/fft_amplitude.png', dpi=300)
            plt.close()

            param.append(float(subdirectory.split("=")[-1]))
        np.savetxt(directory + '/param.txt', param, delimiter=',')
        np.savetxt(directory + '/freq.txt', freq, delimiter=',')
        np.savetxt(directory + '/ampd.txt', ampd, delimiter=',')
        np.savetxt(directory + '/power.txt', av_power, delimiter=',')

        np.savetxt(process_directory + '/param_' + str(sigma_i) + '.txt', param, delimiter=',')
        np.savetxt(process_directory + '/freq_' + str(sigma_i) + '.txt', freq, delimiter=',')
        np.savetxt(process_directory + '/ampd_' + str(sigma_i) + '.txt', ampd, delimiter=',')
        np.savetxt(process_directory + '/power_' + str(sigma_i) + '.txt', av_power, delimiter=',')

        plt.scatter(param, freq, c="k", zorder=10)
        plt.grid(alpha=0.5, zorder=1)
        plt.ylim([0, 1.1 * np.amax(freq)])
        plt.xlabel("$\gamma_0$", fontsize=20)
        plt.ylabel("$f$", fontsize=20)
        plt.tight_layout()
        plt.savefig(directory + '/freq_oscillations.png', dpi=300)
        plt.close()

        plt.scatter(param, ampd, c="k", zorder=10)
        plt.grid(alpha=0.5, zorder=1)
        plt.ylim([0, 1.1 * np.amax(ampd)])
        plt.xlabel("$\gamma_0$", fontsize=20)
        plt.ylabel("$\delta A_R$", fontsize=20)
        plt.tight_layout()
        plt.savefig(directory + '/ampd_oscillations.png', dpi=300)
        plt.close()

        plt.scatter(param, av_power, c="k", zorder=10)
        plt.grid(alpha=0.5, zorder=1)
        plt.ylim([0, 1.1 * np.amax(av_power)])
        plt.xlabel("$\gamma_0$", fontsize=20)
        plt.ylabel("$\\textrm{Average Power}$", fontsize=20)
        plt.tight_layout()
        plt.savefig(directory + '/power_oscillations.png', dpi=300)
        plt.close()