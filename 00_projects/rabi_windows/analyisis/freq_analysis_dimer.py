from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    frequencies = []
    powers = []
    modules = []
    disc = 'D:/'

    initial_dir_data = str(disc) + 'Users/mnustes_science/PT_fluids/mnustes_science/simulation_data'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    distances = np.loadtxt(directory + '/analysis/dists.txt', delimiter=',')
    T = np.loadtxt(directory + '/analysis/t_grid.txt', delimiter=',')
    UR_Rs = np.loadtxt(directory + '/analysis/U1_Rs.txt', delimiter=',')
    UR_Is = np.loadtxt(directory + '/analysis/U1_Is.txt', delimiter=',')
    UL_Rs = np.loadtxt(directory + '/analysis/U2_Rs.txt', delimiter=',')
    UL_Is = np.loadtxt(directory + '/analysis/U2_Is.txt', delimiter=',')
    save_directory = directory + '/analysis'

    Nt = len(T)
    Nt_initial = int(0.2 * Nt)
    T = T[Nt_initial:] - T[Nt_initial]
    for i in range(len(distances)):
        print("###### d = " + str(distances[i]) + " ######")
        dt = T[1] - T[0]
        CCF = np.correlate(UR_Rs[i, Nt_initial:], UR_Rs[i, Nt_initial:], "full")
        tau = np.arange(-T[-1], T[-1], dt)
        Ntau = len(tau)
        dtau = tau[1] - tau[0]

        CCF_max, tau_max, CCF_I = max_finder(CCF[:-1], tau, Ntau, dtau)
        tau_R = []
        maxval = np.amax(CCF)
        power = np.sum(np.abs(UR_Rs[i, Nt_initial:]) ** 2 + np.abs(UR_Is[i, Nt_initial:]) ** 2) / (T[-1] - T[0])
        UR = UR_Rs[i, Nt_initial:] + 1j * UR_Is[i, Nt_initial:]
        UL = UL_Rs[i, Nt_initial:] + 1j * UL_Is[i, Nt_initial:]
        phi_1 = (UL + 1j * UR) / 2
        phi_2 = (UL - 1j * UR) / 2
        R = np.mean(np.abs(phi_1))
        R_std = np.std(np.abs(phi_1))
        P = np.mean(np.abs(phi_2))
        P_std = np.std(np.abs(phi_2))
        for j in range(len(tau_max)):
            if CCF_max[j] > 0.25 * maxval:
                tau_R.append(tau_max[j])
        if len(tau_R) == 1:
            freq = 0
            freq_std = 0
        else:
            freq = 1 / np.mean(np.diff(tau_R))
            freq_std = 1 * np.abs(np.std(np.diff(tau_R)) / np.mean(np.diff(tau_R)) ** 2)
        frequencies.append([freq, freq_std])
        powers.append(power)
        modules.append([R, R_std, P, P_std])
    frequencies = np.array(frequencies)
    powers = np.array(powers)
    modules = np.array(modules)
    np.savetxt(save_directory + '/freqs.txt', frequencies, delimiter=',')
    np.savetxt(save_directory + '/powers.txt', powers, delimiter=',')
    np.savetxt(save_directory + '/modules.txt', modules, delimiter=',')
    plt.plot(distances, frequencies)
    plt.show()