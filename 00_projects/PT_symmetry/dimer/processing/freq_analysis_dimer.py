from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    frequencies = []
    disc = "C:"
    nu = "0.020"
    gamma = "0.185"
    sigma = "6.000"
    alpha = "6.524"

    distances = np.loadtxt(disc + '/mnustes_science/simulation_data/FD/PT_dimer/alpha=' + alpha + '/beta=1.000/mu=0.100/nu=' + nu +'/sigma=' + sigma + '/gamma=' + gamma + '/analysis/dists.txt', delimiter=',')
    T = np.loadtxt(disc + '/mnustes_science/simulation_data/FD/PT_dimer/alpha=' + alpha + '/beta=1.000/mu=0.100/nu=' + nu +'/sigma=' + sigma + '/gamma=' + gamma + '/analysis/t_grid.txt', delimiter=',')
    UR_Rs = np.loadtxt(disc + '/mnustes_science/simulation_data/FD/PT_dimer/alpha=' + alpha + '/beta=1.000/mu=0.100/nu=' + nu +'/sigma=' + sigma + '/gamma=' + gamma + '/analysis/UR_Rs.txt', delimiter=',')
    UR_Is = np.loadtxt(disc + '/mnustes_science/simulation_data/FD/PT_dimer/alpha=' + alpha + '/beta=1.000/mu=0.100/nu=' + nu +'/sigma=' + sigma + '/gamma=' + gamma + '/analysis/UR_Is.txt', delimiter=',')
    UL_Rs = np.loadtxt(disc + '/mnustes_science/simulation_data/FD/PT_dimer/alpha=' + alpha + '/beta=1.000/mu=0.100/nu=' + nu +'/sigma=' + sigma + '/gamma=' + gamma + '/analysis/UL_Rs.txt', delimiter=',')
    UL_Is = np.loadtxt(disc + '/mnustes_science/simulation_data/FD/PT_dimer/alpha=' + alpha + '/beta=1.000/mu=0.100/nu=' + nu +'/sigma=' + sigma + '/gamma=' + gamma + '/analysis/UL_Is.txt', delimiter=',')
    save_directory = disc + '/mnustes_science/simulation_data/FD/PT_dimer/alpha=' + alpha + '/beta=1.000/mu=0.100/nu=' + nu +'/sigma=' + sigma + '/gamma=' + gamma + '/analysis'

    for i in range(len(distances)):
        print(distances[i])
        dt = T[1] - T[0]
        CCF = np.correlate(UR_Rs[i], UR_Rs[i], "full")
        tau = np.arange(-T[-1], T[-1], dt)

        Ntau = len(tau)
        dtau = tau[1] - tau[0]

        CCF_max, tau_max, CCF_I = max_finder(CCF[:-1], tau, Ntau, dtau)
        #plt.plot(tau, CCF[:-1])
        #plt.scatter(tau_max, CCF_max)
        #plt.show()
        #plt.close()
        tau_R = []
        maxval = np.amax(CCF)
        for j in range(len(tau_max)):
            if CCF_max[j] > 0.25 * maxval:
                tau_R.append(tau_max[j])
        if len(tau_R) == 1:
            freq = 0
            freq_std = 0
        else:
            freq = 2 * np.pi * 1000 / np.mean(np.diff(tau_R))
            freq_std = 2 * np.pi * 1000 * np.abs(np.std(np.diff(tau_R)) / np.mean(np.diff(tau_R)) ** 2)
        frequencies.append([freq, freq_std])
    frequencies = np.array(frequencies)
    np.savetxt(save_directory + '/freqs.txt', frequencies, delimiter=',')