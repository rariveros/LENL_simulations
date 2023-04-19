from directories_lyap import *

if __name__ == '__main__':
    lyap_mean_01 = np.loadtxt("E:\mnustes_science\simulation_data\FD\PDNLS_chaos\simulations\lug_lef\E=06\lyap.txt", delimiter=',')
    #lyap_mean_02 = np.loadtxt("E:\mnustes_science\simulation_data\FD\PDNLS_chaos\simulations\lug_lef\E=10\lyap.txt", delimiter=',')
    plt.plot(np.flip(np.sort(lyap_mean_01)), label="$E_i = 5$", linewidth=0.5)
    #plt.plot(np.flip(np.sort(lyap_mean_02)), label="$E_i = 10$", linewidth=0.5)
    plt.title("$\\textrm{Lyapunov Spectrum}$", size='20')
    plt.legend(fontsize=8)
    plt.xlabel('$\\textrm{Initial Conditions}$', size='20')
    plt.xticks(fontsize=15)
    plt.xlim([0, 400])
    plt.ylabel('$\lambda$', size='20')
    plt.yticks(fontsize=15)
    plt.ylim([-1, 0.3])
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig('E:\mnustes_science\simulation_data\FD\PDNLS_chaos\simulations\lug_lef\lyap_espectrums_zoom.png', dpi=300)
    plt.close()