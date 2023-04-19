from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    disco = 'E:/'
    initial_dir_data = str(disco) + "mnustes_science/simulation_data/FD/PDNLS_chaos/branching_dynamics"
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    directories = os.listdir(working_directory)
    for directory_i in directories:
        directory = working_directory + "/" + directory_i
        print(directory_i)

        Z_r = np.loadtxt(directory + '/field_real_filt.txt', delimiter=',')
        Z_i = np.loadtxt(directory + '/field_img_filt.txt', delimiter=',')
        X = np.loadtxt(directory + '/X.txt', delimiter=',')
        T = np.loadtxt(directory + '/T.txt', delimiter=',')

        Nx = len(X)
        Nt = len(T)
        #ti = 0
        #tf = -1

        #Z_r = Z_r[ti:tf, :]
        #T = T[ti:tf]
        Nt = len(T)
        dt = T[1] - T[0]
        dx = X[1] - X[0]

        Nx = len(X)
        N = Nx
        Z_r_fft = []
        for i in range(Nt):
            yf_i = fft(Z_r[i, :])
            yf_i = 2.0 / N * np.abs(yf_i[0:N // 2])
            xf_i = fftfreq(N, dx)[:N // 2]
            Z_r_fft.append(yf_i.tolist())
        K = xf_i
        Z_r_fft = np.array(Z_r_fft)

        T_fft = T

        np.savetxt(directory + '/field_real_FFT.txt', Z_r_fft, delimiter=',')
        np.savetxt(directory + '/K.txt', K, delimiter=',')
        np.savetxt(directory + '/T.txt', T_fft, delimiter=',')

        fig, ax = plt.subplots()

        # legend
        pcm = plt.pcolormesh(xf_i, T, Z_r_fft, cmap="jet", shading='auto')
        cbar = plt.colorbar(pcm)
        cbar.set_label('$\hat{A}_R(t, k)$', rotation=0, size=20, labelpad=-27, y=1.11)

        # put the major ticks at the middle of each cell
        plt.xlim([0, 0.6])
        plt.xlabel('$k$', size='20')
        plt.xticks(fontsize=15)

        plt.ylabel('$T$', size='20')
        plt.yticks(fontsize=15)
        #plt.ylim([0, 0.5])

        plt.grid(linestyle='--', alpha=0.2, color='k')
        cbar.ax.tick_params(labelsize=15)

        # labels
        plt.tight_layout()
        plt.savefig(directory + '/spatial_freq.png', dpi=300)
        plt.close()
        print("--------------------------------------")