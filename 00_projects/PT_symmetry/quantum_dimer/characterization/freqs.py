from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    frequencies = []
    disc = "D:"
    initial_dir_data = str(disc) + 'mnustes_science/simulation_data/FD'
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    directories = [name for name in os.listdir(working_directory) if os.path.isdir(os.path.join(working_directory, name))]
    FREQS = []
    OMEGAS = []
    K = []
    for directory_01 in directories:
        dir_01 = working_directory + "/" + directory_01
        directories_01 = [name for name in os.listdir(dir_01) if os.path.isdir(os.path.join(dir_01, name))]
        frequencies = []
        omegas = []
        ks = []
        for directory_02 in directories_01:
            dir_02 = working_directory + "/" + directory_01 + "/" + directory_02
            params = np.loadtxt(dir_02 + '/parameters.txt', delimiter=',')
            T = np.loadtxt(dir_02 + '/T.txt', delimiter=',')
            U1 = np.loadtxt(dir_02 + '/U1.txt', delimiter=',', dtype=np.complex128)
            V1 = np.loadtxt(dir_02 + '/V1.txt', delimiter=',', dtype=np.complex128)
            save_directory = working_directory
            Nt = len(T)
            t0 = int(0.2 * Nt)

            T = T[t0:] - T[t0]
            U1 = U1[t0:]
            V1 = V1[t0:]
            Nt = len(T)
            dt = T[1] - T[0]
            CCF = np.correlate(np.abs(np.imag(U1)), np.abs(np.imag(U1)  ), "full")
            tau = np.arange(-T[-1] + dt, T[-1] - dt, dt)
            Ntau = len(tau)
            dtau = tau[1] - tau[0]

            CCF_max, tau_max, CCF_I = max_finder(CCF[:-1], tau, Ntau, dtau)
            tau_R = []
            maxval = np.amax(CCF)
            for j in range(len(tau_max)):
                if CCF_max[j] > 0.5 * maxval:
                    tau_R.append(tau_max[j])
            if np.std(np.imag(U1)) < 0.01:
                freq = 0
                freq_std = 0
            else:
                freq = 1 / np.mean(np.diff(tau_R))
                freq_std = 1 * np.abs(np.std(np.diff(tau_R)) / np.mean(np.diff(tau_R)) ** 2)

            # [Delta, gamma, Omega, k, g]
            omega = params[2]
            k = params[3]
            print("## Omega = " + str(omega) + " #### k = " + str(k))
            print(np.std(np.imag(U1)))
            #frequencies.append([freq, freq_std])
            frequencies.append(freq)
            ks.append(k)
            omegas.append(omega)
        FREQS.append(frequencies)
        K.append(ks)
        OMEGAS.append(omegas)
    FREQS = np.array(FREQS)
    K = np.array(K)
    OMEGAS = np.array(OMEGAS)

    np.savetxt(save_directory + '/frequencies.txt', FREQS, delimiter=',')
    np.savetxt(save_directory + '/ks.txt', K, delimiter=',')
    np.savetxt(save_directory + '/omegas.txt', OMEGAS, delimiter=',')

