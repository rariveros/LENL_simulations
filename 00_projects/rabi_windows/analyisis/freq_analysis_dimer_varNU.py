from functions import *
from back_process import *
from time_integrators import *

def gaussian(x, A, x_0, sigma):
    return A * np.exp(- (x - x_0) ** 2 / (2 * sigma ** 2))


def resample(old_arrays, old_x, N_resample):
    new_arrays = []
    for i in range(len(old_arrays)):
        [new_array_i, new_x] = signal.resample(old_arrays[i], int(len(old_x) * (N_resample)), old_x)
        new_arrays.append(new_array_i)
    return new_arrays, new_x


if __name__ == '__main__':

    eq = 'PT_dimer'
    gamma_str = "0.280"
    sigme_str = "3.000"
    disco = 'D:/'
    initial_dir_data = str(disco) + 'Users/mnustes_science/PT_fluids/mnustes_science/simulation_data'
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    directories = [name for name in os.listdir(working_directory) if os.path.isdir(os.path.join(working_directory, name))]
    for directory in directories:
        print("##### " + str(directory) + " #####")
        frequencies = []
        powers = []
        modules = []

        distances = np.loadtxt(working_directory + '/' + directory + '/sigma=' + sigme_str +'/gamma=' + gamma_str + '/analysis/dists.txt', delimiter=',')
        T = np.loadtxt(working_directory + '/' + directory + '/sigma=' + sigme_str +'/gamma=' + gamma_str + '/analysis/t_grid.txt', delimiter=',')
        UR_Rs = np.loadtxt(working_directory + '/' + directory + '/sigma=' + sigme_str +'/gamma=' + gamma_str + '/analysis/U1_Rs.txt', delimiter=',')
        UR_Is = np.loadtxt(working_directory + '/' + directory + '/sigma=' + sigme_str +'/gamma=' + gamma_str + '/analysis/U1_Is.txt', delimiter=',')
        UL_Rs = np.loadtxt(working_directory + '/' + directory + '/sigma=' + sigme_str +'/gamma=' + gamma_str + '/analysis/U2_Rs.txt', delimiter=',')
        UL_Is = np.loadtxt(working_directory + '/' + directory + '/sigma=' + sigme_str +'/gamma=' + gamma_str + '/analysis/U2_Is.txt', delimiter=',')
        save_directory = working_directory + '/' + directory + '/sigma=' + sigme_str +'/gamma=' + gamma_str +  '/analysis'

        Nt = len(T)
        Nt_initial = int(0.2 * Nt)
        T = T[Nt_initial:] - T[Nt_initial]
        for i in range(len(distances)):
            print("## d = " + str(distances[i]) + " ##")
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