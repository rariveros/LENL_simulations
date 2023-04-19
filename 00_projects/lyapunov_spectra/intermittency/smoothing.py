from functions import *
from back_process import *
from time_integrators import *
from directories_lyap import *

if __name__ == '__main__':
    disco = 'E:/'
    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FD' + main_directory
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    directories = os.listdir(working_directory)
    for directory_i in directories:
        directory = working_directory + "/" + directory_i
        Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
        Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
        X = np.loadtxt(directory + '/X.txt', delimiter=',')
        T = np.loadtxt(directory + '/T.txt', delimiter=',')

        Z_r_filtered = filtro_superficie(Z_r, 200, "Y")
        Z_i_filtered = filtro_superficie(Z_r, 200, "Y")

        Z_r_fine = Z_r - Z_r_filtered
        Z_i_fine = Z_i - Z_i_filtered
        [alpha, beta, gamma_0, mu_0, nu] = np.loadtxt(directory + '/parameters.txt', delimiter=',')
        sigma = directory_i.split("=")
        sigma = float(sigma[-1])
        print("sigma = " + str(sigma))
        subfile = pdnlS_name([alpha, beta, mu_0, nu, sigma, gamma_0], "ABMNGS")
        save_directory = initial_dir_data + subdir_branching + subfile
        if not os.path.exists(save_directory):
            os.makedirs(save_directory)
        np.savetxt(save_directory + '/field_real_fine.txt', Z_r_fine, delimiter=',')
        np.savetxt(save_directory + '/field_img_fine.txt', Z_i_fine, delimiter=',')
        np.savetxt(save_directory + '/field_real_filt.txt', Z_r_filtered, delimiter=',')
        np.savetxt(save_directory + '/field_img_filt.txt', Z_i_filtered, delimiter=',')
        np.savetxt(save_directory + '/X.txt', X, delimiter=',')
        np.savetxt(save_directory + '/T.txt', T, delimiter=',')

        ####

        fig, ax = plt.subplots()

        # legend
        pcm = plt.pcolormesh(X, T, Z_r_fine, cmap=parula_map, vmin=np.amin(np.real(Z_r_fine)), vmax=np.amax(np.real(Z_r_fine)), shading='auto')
        cbar = plt.colorbar(pcm)
        cbar.set_label('$A_R(x, t)$', rotation=0, size=20, labelpad=-27, y=1.11)

        # put the major ticks at the middle of each cell
        plt.xlim([X[0], X[-1]])
        plt.xlabel('$x$', size='15')
        plt.xticks(fontsize=20)

        #plt.ylim([20, 80])
        plt.ylabel('$t$', size='15')
        plt.yticks(fontsize=20)
        #plt.ylim([0, 0.5])

        plt.grid(linestyle='--', alpha=0.2, color='k')
        cbar.ax.tick_params(labelsize=15)

        # labels
        plt.tight_layout()
        plt.savefig(save_directory + '/field_real_fine.png', dpi=300)
        plt.close()

        ####

        fig, ax = plt.subplots()

        # legend
        pcm = plt.pcolormesh(X, T, Z_r_filtered, cmap=parula_map, vmin=np.amin(np.real(Z_r_filtered)), vmax=np.amax(np.real(Z_r_filtered)), shading='auto')
        cbar = plt.colorbar(pcm)
        cbar.set_label('$A_R(x, t)$', rotation=0, size=20, labelpad=-27, y=1.11)

        # put the major ticks at the middle of each cell
        plt.xlim([X[0], X[-1]])
        plt.xlabel('$x$', size='15')
        plt.xticks(fontsize=20)

        #plt.ylim([20, 80])
        plt.ylabel('$t$', size='15')
        plt.yticks(fontsize=20)
        #plt.ylim([0, 0.5])

        plt.grid(linestyle='--', alpha=0.2, color='k')
        cbar.ax.tick_params(labelsize=15)

        # labels
        plt.tight_layout()
        plt.savefig(save_directory + '/field_real_filt.png', dpi=300)
        plt.close()
        print("----------------------------------")