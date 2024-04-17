from functions import *
from back_process import *
from time_integrators import *
from scipy import stats


if __name__ == '__main__':

    # Definiendo parámetros
    directory = "C:/mnustes_science/simulation_data/FD/localized_chaos/EE_02"
    PDF_EE = np.loadtxt(directory + "/PDF_EEs.txt", delimiter=',')
    dH_EE = np.loadtxt(directory + "/dH_EEs.txt", delimiter=',')
    gammas = np.loadtxt(directory + "/gammas.txt", delimiter=',')
    p_EE = np.loadtxt(directory + "/p_EE.txt", delimiter=',')
    n = 0
    for i in range(len(dH_EE)):
        mean_i = np.sum(dH_EE[i] * np.exp(PDF_EE[i])) / np.sum(np.exp(PDF_EE[i]))
        gamma_i = gammas[i+1]
        spectrum = PDF_EE[i]
        spectrum = filtro_array(3, spectrum)
        plt.plot(dH_EE[i] / mean_i, np.log(np.exp(spectrum) / np.sum(np.exp(PDF_EE[i]))), color=(1 - (n / 13), 0, n / 13), label="$\gamma_0 = " + f"{gamma_i:.3f}" + "$")
        n = n + 1
    plt.xlabel('$h/\langle h \\rangle$', size='20')
    plt.xticks(fontsize=15)
    plt.xlim(0, 5)
    plt.yticks(fontsize=15)
    plt.ylabel('$\\textrm{ln(PDF)}$', size='20')
    plt.ylim([-12, 0])
    plt.grid(alpha=0.2)
    plt.legend()
    plt.tight_layout()
    plt.savefig('EE.png', dpi=300)
    plt.show()
    plt.close()


