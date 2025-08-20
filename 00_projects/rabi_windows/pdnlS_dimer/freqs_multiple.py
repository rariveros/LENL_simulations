from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    disc = "D:"
    initial_dir_data = str(disc) + 'mnustes_science/simulation_data/FD'
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    directories = [name for name in os.listdir(working_directory) if os.path.isdir(os.path.join(working_directory, name))]
    for directory in directories:
        dir = working_directory + "/" + directory + "/mu=0.1000/gamma=0.2000"
        subdirectories = [name for name in os.listdir(dir) if os.path.isdir(os.path.join(dir, name))]
        FREQS = []
        K = []
        nus = []
        modules = []
        for directory_01 in subdirectories:
            dir_01 = dir + "/" + directory_01
            params = np.loadtxt(dir_01 + '/parameters.txt', delimiter=',')
            T = np.loadtxt(dir_01 + '/T.txt', delimiter=',')
            U = np.loadtxt(dir_01 + '/U.txt', delimiter=',', dtype=np.complex128)
            V = np.loadtxt(dir_01 + '/V.txt', delimiter=',', dtype=np.complex128)
            save_directory = dir
            Nt = len(T)
            t0 = int(0.2 * Nt)

            T = T[t0:] - T[t0]
            U = U[t0:]
            V = V[t0:]
            Nt = len(T)
            dt = T[1] - T[0]
            CCF = np.correlate(np.abs(np.imag(U)), np.abs(np.imag(U)), "full")
            tau = np.arange(-T[-1] + dt, T[-1] - dt, dt)
            Ntau = len(tau)
            dtau = tau[1] - tau[0]

            CCF_max, tau_max, CCF_I = max_finder(CCF[:-1], tau, Ntau, dtau)
            tau_R = []
            maxval = np.amax(CCF)
            R = np.mean(np.abs(U))
            R_std = np.std(np.abs(U))
            P = np.mean(np.abs(V))
            P_std = np.std(np.abs(V))
            for j in range(len(tau_max)):
                if CCF_max[j] > 0.25 * maxval:
                    tau_R.append(tau_max[j])
            if np.std(np.abs(np.imag(U))) < 0.001:
                freq = 0
                freq_std = 0
            else:
                freq = 1 / np.mean(np.diff(tau_R))
                freq_std = 1 * np.abs(np.std(np.diff(tau_R)) / np.mean(np.diff(tau_R)) ** 2)

            [nu, mu, gamma, k, g] = params
            print("#### k = " + str(k) + " ####")
            FREQS.append(freq)
            modules.append([R, R_std, P, P_std])
            K.append(k)
            nus.append(nu)
        FREQS = np.array(FREQS)
        K = np.array(K)
        nus = np.array(nus)
        modules = np.array(modules)
        np.savetxt(save_directory + '/frequencies.txt', np.array(FREQS), delimiter=',')
        np.savetxt(save_directory + '/Ks.txt', K, delimiter=',')
        np.savetxt(save_directory + '/nus.txt', nus, delimiter=',')
        np.savetxt(save_directory + '/modules.txt', modules, delimiter=',')

        plt.scatter(K, FREQS, c="k")
        plt.xlabel("$\kappa$")
        plt.ylabel("$\Omega$")
        plt.savefig(save_directory + "/freqs.png", dpi=300)
        plt.close()
