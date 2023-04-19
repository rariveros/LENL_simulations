from functions import *
from back_process import *
from time_integrators import *
from directories_lyap import *

if __name__ == '__main__':
    disco = 'E:/'
    initial_dir_data = str(disco) + "mnustes_science/simulation_data/FD/PDNLS_chaos"
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    directories = os.listdir(working_directory)
    for directory_i in directories:
        directory = working_directory + "/" + directory_i
        print(directory_i)

        sigma = float(directory.split("=")[-1])
        print(sigma)
        [alpha, beta, gamma_0, mu_0, nu] = np.loadtxt(directory + '/parameters.txt', delimiter=',')

        Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
        Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
        X = np.loadtxt(directory + '/X.txt', delimiter=',')
        T = np.loadtxt(directory + '/T.txt', delimiter=',')

        Nx = len(X)
        Nt = len(T)
        ti = 0
        tf = -1

        Z_r = Z_r[ti:tf, :]
        Z_i = Z_i[ti:tf, :]
        T = T[ti:tf]
        Nt = len(T)
        dt = T[1] - T[0]
        dx = X[1] - X[0]

        Nx = len(X)
        N = Nx
        Z_r = np.transpose(Z_r)
        Z_i = np.transpose(Z_i)
        D = sparse_D_neumann(Nx, dx)
        hamiltonian = Der(D, Z_r) ** 2 + Der(D, Z_i) ** 2 - 0.5 * (Z_r ** 2 + Z_i ** 2) ** 2 - nu * (Z_r ** 2 + Z_i ** 2)
        hamiltonian = np.transpose(hamiltonian)
        Hamiltonian = []
        for i in range(Nt):
            int_i = integrate.simpson(hamiltonian[i, :], X)
            Hamiltonian.append(int_i)
        H = np.array(Hamiltonian)

        subfile = pdnlS_name([alpha, beta, mu_0, nu, sigma, gamma_0], "ABMNGS")
        save_directory = initial_dir_data + subdir_branching + subfile
        if not os.path.exists(save_directory):
            os.makedirs(save_directory)
        H_filt = filtro_array(300, H)
        np.savetxt(save_directory + '/hamiltonian_xt.txt', hamiltonian, delimiter=',')
        np.savetxt(save_directory + '/Hamiltonian_t.txt', H, delimiter=',')
        np.savetxt(save_directory + '/Hamiltonian_averaged.txt', H_filt, delimiter=',')


        plt.plot(T, H, zorder=0, alpha=0.5, color="r", label="$H(t)$")
        plt.plot(T, H_filt, zorder=1, color="k", label="$\langle H(t) \\rangle_{t}$")
        plt.legend()
        plt.xlabel('$t$', size='20')
        plt.xticks(fontsize=15)
        plt.xlim([T[0], T[-1]])
        plt.ylabel('$H(t)$', size='20')
        plt.yticks(fontsize=15)
        plt.grid(alpha=0.2)
        plt.tight_layout()
        plt.savefig(save_directory + '/hamiltonian_mean.png', dpi=300)
        plt.close()

        #fig, ax = plt.subplots()

        # legend
        #pcm = plt.pcolormesh(X, T, hamiltonian, cmap="jet", shading='auto')
        #cbar = plt.colorbar(pcm)
        #cbar.set_label('$\mathcal{H}(x, t)$', rotation=0, size=20, labelpad=-27, y=1.11)

        # put the major ticks at the middle of each cell
        #plt.xlabel('$k$', size='20')
        #plt.xticks(fontsize=15)

        #plt.ylabel('$T$', size='20')
        #plt.yticks(fontsize=15)

        #plt.grid(linestyle='--', alpha=0.2, color='k')
        #cbar.ax.tick_params(labelsize=15)

        # labels
        #plt.tight_layout()
        #plt.savefig(save_directory + '/hamiltonian.png', dpi=300)
        #plt.close()
        print("----------------------------------------------")