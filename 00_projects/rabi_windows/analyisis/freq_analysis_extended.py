from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    disco = 'D:/'
    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FDa'
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    directories = [name for name in os.listdir(working_directory) if os.path.isdir(os.path.join(working_directory, name))]
    ### E:\mnustes_science\simulation_data\FD\rabi_windows\dimensionless\var_dist\alpha=1.000\beta=1.000\mu=0.1000\nu=0.3200\sigma=3.000\gamma=0.2800
    freqs = []
    powers = []
    dists = []
    nus = []
    for directory in directories:
        print("#############   " + directory + "   #############")
        save_directory = working_directory + "/"

        Z_r = np.loadtxt(working_directory + "/" + directory + '/field_real.txt', delimiter=',')
        Z_i = np.loadtxt(working_directory + "/" + directory + '/field_img.txt', delimiter=',')
        X = np.loadtxt(working_directory + "/" + directory + '/X.txt', delimiter=',')
        T = np.loadtxt(working_directory + "/" + directory + '/T.txt', delimiter=',')
        params = np.loadtxt(working_directory + "/" + directory + '/parameters.txt', delimiter=',')
        #[alpha, beta, gamma_0, dist,  mu, nu, sigma]
        dist = params[3]
        nu = params[6]
        Nx = len(X)
        Nt = len(T)

        Nt_init = int(0.3 * Nt)
        Z_r = Z_r[Nt_init:]
        Z_i = Z_i[Nt_init:]
        T = T[Nt_init:] - T[Nt_init]

        ZR_L = Z_r[:, 0:int(Nx / 2)]
        ZR_R = Z_r[:, int(Nx / 2) + 1:]
        ZI_L = Z_i[:, 0:int(Nx / 2)]
        ZI_R = Z_i[:, int(Nx / 2) + 1:]
        X_L = X[0:int(Nx / 2)]
        X_R = X[int(Nx / 2) + 1:]
        UR_R = integrate.simpson(ZR_R, X_R)
        UR_L = integrate.simpson(ZR_L, X_L)
        UI_R = integrate.simpson(ZI_R, X_R)
        UI_L = integrate.simpson(ZI_L, X_L)

        dt = T[1] - T[0]
        CCF = np.correlate(UR_R, UR_R, "full")
        tau = np.arange(-T[-1], T[-1], dt)

        Ntau = len(tau)
        dtau = tau[1] - tau[0]
        CCF_max, tau_max, CCF_I = max_finder(CCF[:-1], tau, Ntau, dtau)
        tau_R = []
        maxval = np.amax(CCF)
        power = np.sum(np.abs(UR_R) ** 2 + np.abs(UI_R) ** 2) / (T[-1] - T[0])
        power_std = np.std(np.abs(UR_R))
        for j in range(len(tau_max)):
            if CCF_max[j] > 0.25 * maxval:
                tau_R.append(tau_max[j])
        if power_std < 0.001 and power < 0.001:
            freq = 0
            freq_std = 0
        elif len(tau_R) == 1:
            freq = 0
            freq_std = 0
        else:
            freq = 1 / np.mean(np.diff(tau_R))
            freq_std = 1 * np.abs(np.std(np.diff(tau_R)) / np.mean(np.diff(tau_R)) ** 2)
        freqs.append([freq, freq_std])
        powers.append([power, power_std])
        dists.append(dist)
        nus.append(nu)
    freqs = np.array(freqs)
    powers = np.array(powers)
    dists = np.array(dists)
    nus = np.array(nus)

    # Get sorting indices from dists
    sort_idx = np.argsort(dists)

    # Apply sorting to all arrays
    freqs = freqs[sort_idx]
    powers = powers[sort_idx]
    dists = dists[sort_idx]
    nus = nus[sort_idx]

    np.savetxt(save_directory + '/freqs.txt', freqs, delimiter=',')
    np.savetxt(save_directory + '/powers.txt', powers, delimiter=',')
    np.savetxt(save_directory + '/dists.txt', dists, delimiter=',')
    np.savetxt(save_directory + '/nus.txt', nus, delimiter=',')