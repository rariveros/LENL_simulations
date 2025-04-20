from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    disc = "aD:"
    initial_dir_data = str(disc) + 'mnustes_science/simulation_data/FD'
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    directories = [name for name in os.listdir(working_directory) if os.path.isdir(os.path.join(working_directory, name))]
    FREQS = []
    OMEGAS = []
    Rs = []
    Ss = []
    THETA = []
    PHI = []
    for directory_01 in directories:
        dir_01 = working_directory + "/" + directory_01
        params = np.loadtxt(dir_01 + '/parameters.txt', delimiter=',')
        T = np.loadtxt(dir_01 + '/T.txt', delimiter=',')
        U1 = np.loadtxt(dir_01 + '/U1.txt', delimiter=',', dtype=np.complex128)
        V2 = np.loadtxt(dir_01 + '/V2.txt', delimiter=',', dtype=np.complex128)
        save_directory = working_directory
        Nt = len(T)
        t0 = int(0.5 * Nt)

        T = T[t0:] - T[t0]
        U1 = U1[t0:]
        V1 = V2[t0:]
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
            if CCF_max[j] > 0.7 * maxval:
                tau_R.append(tau_max[j])
        if np.std(np.imag(U1)) < 0.05:
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
        FREQS.append(freq)
        OMEGAS.append(omega)
        Rs.append(np.amax(np.abs(U1)))
        Ss.append(np.amax(np.abs(V1)))
        theta = np.arctan(np.amax(np.abs(V1)) * np.cos(np.angle(U1[-1])) / (np.amax(np.abs(U1)) * np.sin(np.angle(V1[-1]))))
        phi = np.arctan(-(np.amax(np.abs(U1)) * np.sin(np.angle(U1[-1]))) / (np.amax(np.abs(V1)) * np.cos(np.angle(V1[-1]))))
        THETA.append(theta)
        PHI.append(phi)
        plt.plot(np.abs(U1), np.abs(V1), color="k")
    plt.show()
    FREQS = np.array(FREQS)
    OMEGAS = np.array(OMEGAS)
    Rs = np.array(Rs)
    Ss = np.array(Ss)
    np.savetxt(save_directory + '/frequencies.txt', np.array(FREQS), delimiter=',')
    np.savetxt(save_directory + '/omegas.txt', OMEGAS, delimiter=',')

    plt.scatter(OMEGAS, FREQS, c="k")
    plt.show()
    plt.close()

    nu = 0.0 + 1j * 0
    mu = 0.1

    R0 = P0 = np.sqrt(-nu + np.sqrt(OMEGAS ** 2 - mu ** 2))
    R1 = (R0 / OMEGAS) * ((np.sqrt(OMEGAS ** 2 - mu ** 2) + mu) / (- 2 * nu + 4 * np.sqrt(OMEGAS ** 2 - mu ** 2)))
    P1 = -(P0 / OMEGAS) * ((np.sqrt(OMEGAS ** 2 - mu ** 2) + mu) / (- 2 * nu + 4 * np.sqrt(OMEGAS ** 2 - mu ** 2)))
    R = R0 + 0.5 * k * R1  # 0.336 #
    P = P0 + 0.5 * k * P1  # 0.23 #

    #R = np.sqrt(-nu + np.sqrt(OMEGAS ** 2 - (mu - k) ** 2))
    #P = np.sqrt(-nu + np.sqrt(OMEGAS ** 2 - (mu + k) ** 2))

    plt.scatter(OMEGAS, Rs, c="b")
    plt.scatter(OMEGAS, Ss, c="r")
    plt.plot(OMEGAS, R, color="b")
    plt.plot(OMEGAS, P, color="r")
    plt.plot(OMEGAS, R0, color="k")
    plt.show()
    plt.close()



    plt.scatter(OMEGAS, THETA, color='b')
    plt.scatter(OMEGAS, PHI, color='r')
    plt.show()

