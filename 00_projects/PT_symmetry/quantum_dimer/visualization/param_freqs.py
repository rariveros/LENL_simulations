from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    directory = "D:/mnustes_science/simulation_data/FD/ladder_operators/test/Delta=0.1000/gamma=0.1000"
    Freqs = np.loadtxt(directory + '/frequencies.txt', delimiter=',')
    Ks = np.loadtxt(directory + '/ks.txt', delimiter=',')
    Omegas = np.loadtxt(directory + '/omegas.txt', delimiter=',')

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    print(Omegas[:, 0])
    print(Ks[0, :])
    for i in range(0, len(Omegas[:, 0])):
        for j in range(0, len(Ks[0, :])):
            if Omegas[i, j] < 0.055: ## DAMPED
                color = "#FFE14C"
            elif Omegas[i, j] > 0.055 and Freqs[i, j] < 0.001: ## STATIONARY
                color = "#d5a3ff"
            elif Omegas[i, j] > 0.055 and Freqs[i, j] > 0.001: ## RABI
                color = "#bef782"
            ax.scatter(Ks[i, j], Omegas[i, j], Freqs[i, j], c=color, zorder=4, alpha=1.0, edgecolors="k")
        plt.plot(Ks[i, :], Omegas[i, :], Freqs[i, :], color="k", zorder=2)
    ax.view_init(30, -167, 0)
    ax.set_xlim(-0.02, 0.5)
    ax.set_ylim(-0.005, 0.2)
    ax.set_zlim(0, 0.17)
    ax.set_xlabel("$\kappa$", fontsize=20)
    ax.set_ylabel("$\Omega$", fontsize=20)
    ax.set_zlabel("$\\textrm{Frequency}\\times 10^{-3}$", fontsize=20)
    ax.set_xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5], ["$0.0$", "$0.1$", "$0.2$", "$0.3$", "$0.4$", "$0.5$"], fontsize=15)
    ax.set_yticks([0.0, 0.05, 0.1, 0.150, 0.2], ["$0.00$", "$0.05$", "$0.10$", "$0.15$", "$0.20$"], fontsize=15)
    ax.set_zticks([0.0, 0.05, 0.1, 0.150], ["$0.00$", "$0.05$", "$0.10$", "$0.15$"], fontsize=15)
    plt.savefig("freq_space_01.png", dpi=300)