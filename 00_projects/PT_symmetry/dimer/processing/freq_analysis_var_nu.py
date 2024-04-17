from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    disc = "C:"
    nus = ["0.020", "0.040", "0.060", "0.080", "0.100", "0.120", "0.140", "0.160", "0.180", "0.200", "0.240", "0.260", "0.280", "0.300", "0.320"]
    gamma = "0.280"
    sigma = "3.000"
    alpha = "1.000"
    FREQ = []
    NUS = []
    for nu in nus:
        NU = float(nu)
        NUS.append(NU)
        distances = np.loadtxt(disc + '/mnustes_science/simulation_data/FD/PT_dimer/alpha=' + alpha + '/beta=1.000/mu=0.100/nu=' + nu +'/sigma=' + sigma + '/gamma=' + gamma + '/analysis/dists.txt', delimiter=',')
        T = np.loadtxt(disc + '/mnustes_science/simulation_data/FD/PT_dimer/alpha=' + alpha + '/beta=1.000/mu=0.100/nu=' + nu +'/sigma=' + sigma + '/gamma=' + gamma + '/analysis/t_grid.txt', delimiter=',')
        UR_Rs = np.loadtxt(disc + '/mnustes_science/simulation_data/FD/PT_dimer/alpha=' + alpha + '/beta=1.000/mu=0.100/nu=' + nu +'/sigma=' + sigma + '/gamma=' + gamma + '/analysis/UR_Rs.txt', delimiter=',')
        UR_Is = np.loadtxt(disc + '/mnustes_science/simulation_data/FD/PT_dimer/alpha=' + alpha + '/beta=1.000/mu=0.100/nu=' + nu +'/sigma=' + sigma + '/gamma=' + gamma + '/analysis/UR_Is.txt', delimiter=',')
        UL_Rs = np.loadtxt(disc + '/mnustes_science/simulation_data/FD/PT_dimer/alpha=' + alpha + '/beta=1.000/mu=0.100/nu=' + nu +'/sigma=' + sigma + '/gamma=' + gamma + '/analysis/UL_Rs.txt', delimiter=',')
        UL_Is = np.loadtxt(disc + '/mnustes_science/simulation_data/FD/PT_dimer/alpha=' + alpha + '/beta=1.000/mu=0.100/nu=' + nu +'/sigma=' + sigma + '/gamma=' + gamma + '/analysis/UL_Is.txt', delimiter=',')
        save_directory = disc + '/mnustes_science/simulation_data/FD/PT_dimer/alpha=' + alpha + '/beta=1.000/mu=0.100/nu=' + nu +'/sigma=' + sigma + '/gamma=' + gamma + '/analysis'
        frequencies = []
        for i in range(len(distances)):
            print(distances[i])
            dt = T[1] - T[0]
            CCF = np.correlate(UR_Rs[i], UR_Rs[i], "full")
            tau = np.arange(-T[-1], T[-1], dt)

            Ntau = len(tau)
            dtau = tau[1] - tau[0]

            CCF_max, tau_max, CCF_I = max_finder(CCF[:-1], tau, Ntau, dtau)
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
            frequencies.append(freq)
        frequencies = np.array(frequencies)
        FREQ.append(frequencies)
    np.savetxt('FREQS.txt', FREQ, delimiter=',')
    np.savetxt('NUS.txt', NUS, delimiter=',')
    np.savetxt('DIST.txt', distances, delimiter=',')

    FREQ = np.array(FREQ)
    FREQ[np.isnan(FREQ)] = 0
    FREQ = filtro_superficie(FREQ, 4, "YX")

    pcm = plt.pcolormesh(distances, np.array(NUS), np.array(FREQ), cmap="jet", shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    plt.xlabel('$d$', size='20')
    plt.ylabel('$\\nu$', size='20')
    cbar.set_label('$\Omega \\times 10^{-3}$', rotation=0, size=15, labelpad=-27, y=1.1)
    plt.grid(linestyle='--', alpha=0.5)
    plt.savefig("freqs_var_nu.png", dpi=250)
    plt.close()
