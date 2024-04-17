from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    ampd_01 = np.array([7, 7.5, 6.8, 6.6, 6.4, 6.2, 6.0, 5.8, 5.6, 5.4, 5.2, 5.0, 7.7, 7.9, 8.1, 8.3])
    freq_01 = np.array([12.5, 12.4, 12.5, 12.7, 12.6, 12.6, 12.6, 12.8, 12.8, 13.0, 13.0, 13.3, 12.4, 12.4, 12.4, 12.3])
    ampd_02 = np.array([7, 6.8, 6.6, 6.4, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9, 9.2])
    freq_02 = np.array([12.6, 12.7, 12.8, 12.8, 12.6, 12.7, 12.6, 12.5, 12.5, 12.4, 12.4, 12.3, 12.3, 12.3, 12.2])
    d = 20
    alpha, beta, nu_01, gamma_01, f_0 = fluid_pdnls_parameters(freq_01, ampd_01, d)
    alpha, beta, nu_02, gamma_02, f_0 = fluid_pdnls_parameters(freq_02, ampd_02, d)

    alpha, beta, nu_03, gamma_03, f_0 = fluid_pdnls_parameters(10, 15, d)
    plt.scatter(nu_01, gamma_01)
    plt.scatter(nu_02, gamma_02)
    #plt.scatter(nu_03, gamma_03)
    plt.xlabel("$\\nu$", fontsize=20)
    plt.ylabel("$\gamma_0$", fontsize=20)
    plt.xlim([-0.3, 0])
    plt.ylim([0, 0.2])
    plt.grid(alpha=0.2)
    plt.show()

    plt.scatter(freq_01, ampd_01)
    plt.scatter(freq_02, ampd_02)
    plt.xlabel("$f$", fontsize=20)
    plt.ylabel("$a$", fontsize=20)
    #plt.xlim([-0.3, 0])
    #plt.ylim([0, 0.2])
    plt.grid(alpha=0.2)
    plt.show()
