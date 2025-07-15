from functions import *
from back_process import *
from time_integrators import *
from scipy.signal import argrelextrema

if __name__ == '__main__':
    frequencies = []
    disc = "aD:"
    initial_dir_data = str(disc) + 'mnustes_science/simulation_data/FD'
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    directories = [name for name in os.listdir(working_directory) if os.path.isdir(os.path.join(working_directory, name))]

    g = 9790
    l_y = 16
    d = 20
    k_y = np.pi / l_y
    k = k_y
    tau = np.tanh(k * d)
    w_1 = np.sqrt(g * k * tau)

    FREQS = []
    FREQS_STD = []
    AMPDS = []
    DISTS = []
    for directory_01 in directories:
        dir_02 = working_directory + "/" + directory_01
        [alpha, beta, gamma_0, dist, sigma,  mu, nu] = np.loadtxt(dir_02 + '/parameters.txt', delimiter=',')
        T = np.loadtxt(dir_02 + '/T.txt', delimiter=',')
        X = np.loadtxt(dir_02 + '/X.txt', delimiter=',')
        Z_r = np.loadtxt(dir_02 + '/field_real.txt', delimiter=',')
        Z_i = np.loadtxt(dir_02 + '/field_img.txt', delimiter=',')
        Z = Z_r + 1j * Z_i
        Z_mod = np.abs(Z)
        save_directory = working_directory

        T = T / (w_1 / np.pi)
        Nt = len(T)
        t0 = int(0.2 * Nt)
        Nx = len(X)

        T = T[t0:] - T[t0]
        Z_mod = Z_mod[t0:, :]
        Z_mod_L = Z_mod[:, 0:int(Nx / 2)]
        Z_mod_R = Z_mod[:, int(Nx / 2) + 1:]
        X_L = X[0:int(Nx / 2)]
        X_R = X[int(Nx / 2) + 1:]
        U1 = integrate.simpson(Z_mod_L ** 2, X_L)
        V1 = integrate.simpson(Z_mod_R ** 2, X_R)

        Nt = len(T)
        dt = T[1] - T[0]
        CCF = np.correlate(U1, U1, "full")
        tau = np.arange(-T[-1], T[-1], dt)
        Ntau = len(tau)
        dtau = tau[1] - tau[0]
        indices = argrelextrema(CCF, np.greater)[0]
        CCF_max, tau_max, CCF_I = CCF[indices], tau[indices], indices
        #CCF_max, tau_max, CCF_I = max_finder(CCF[:], tau, Ntau, dtau)
        tau_R = []
        maxval = np.amax(CCF)
        for j in range(len(tau_max)):
            if CCF_max[j] > 0.2 * maxval:
                tau_R.append(tau_max[j])
        if np.std(U1) < 0.1:
            freq = 0
            freq_std = 0
        else:
            freq = 1 / np.mean(np.diff(tau_R))
            freq_std = 1 * np.abs(np.std(np.diff(tau_R)) / np.mean(np.diff(tau_R)) ** 2)
        beta = 0.004811649356064012
        FREQS.append(freq)
        FREQS_STD.append(freq_std)
        AMPDS.append(np.amax(Z_mod) / np.sqrt(beta))
        DISTS.append(dist)
    FREQS = np.array(FREQS)
    DISTS = np.array(DISTS)
    np.savetxt(save_directory + '/frequencies.txt', FREQS, delimiter=',')
    np.savetxt(save_directory + '/dists.txt', DISTS, delimiter=',')
    np.savetxt(save_directory + '/freq_std.txt', FREQS_STD, delimiter=',')
    np.savetxt(save_directory + '/ampds.txt', AMPDS, delimiter=',')