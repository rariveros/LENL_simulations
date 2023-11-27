from back_process import *

if __name__ == '__main__':
    sigmas = np.arange(30, 85, 5)
    gammas = np.arange(0.15, 0.21, 0.004)
    for i in range(len(gammas)):
        plt.scatter(gammas[i] * np.ones(len(sigmas)), sigmas, c="k", zorder=10)
    plt.grid(alpha=0.5, zorder=1)
    plt.xlabel("$\gamma_0$", fontsize=25)
    plt.ylabel("$\sigma_i$", fontsize=25)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tight_layout()
    plt.savefig('parameter_space.png', dpi=300)
    plt.close()
