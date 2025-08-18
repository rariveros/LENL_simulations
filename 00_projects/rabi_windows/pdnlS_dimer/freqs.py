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
    FREQS = []
    K = []
    for directory_01 in directories:
        dir_01 = working_directory + "/" + directory_01
        params = np.loadtxt(dir_01 + '/parameters.txt', delimiter=',')
        T = np.loadtxt(dir_01 + '/T.txt', delimiter=',')
        U = np.loadtxt(dir_01 + '/U.txt', delimiter=',', dtype=np.complex128)
        V = np.loadtxt(dir_01 + '/V.txt', delimiter=',', dtype=np.complex128)
        save_directory = working_directory
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
        for j in range(len(tau_max)):
            if CCF_max[j] > 0.7 * maxval:
                tau_R.append(tau_max[j])
        if np.std(np.imag(U)) < 0.05:
            freq = 0
            freq_std = 0
        else:
            freq = 1 / np.mean(np.diff(tau_R))
            freq_std = 1 * np.abs(np.std(np.diff(tau_R)) / np.mean(np.diff(tau_R)) ** 2)

        [nu, mu, gamma, k, g] = params
        print("#### k = " + str(k) + " ####")
        print(np.std(np.imag(U)))
        FREQS.append(freq)
        K.append(k)
    FREQS = np.array(FREQS)
    K = np.array(K)
    np.savetxt(save_directory + '/frequencies.txt', np.array(FREQS), delimiter=',')
    np.savetxt(save_directory + '/Ks.txt', K, delimiter=',')

    plt.scatter(K, FREQS, c="k")
    plt.xlabel("$\kappa$")
    plt.ylabel("$\Omega$")
    plt.savefig(save_directory + "freqs.png", dpi=300)
    plt.close()

