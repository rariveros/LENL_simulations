from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    phase = np.arange(- np.pi, np.pi, 0.01)
    gamma_max = 0.5
    gammas = np.arange(0.1, gamma_max, 0.02)
    OMEGAS_1 = []
    OMEGAS_2 = []
    for gamma_i in gammas:
        gamma_0 = gamma_i
        mu = 0.1
        alpha = 6
        sigma = 90
        amp_i = ((2 * mu / 9) * (gamma_i - mu)) ** 0.25
        A_0 = amp_i * np.exp(1j * phase)
        nu = 0.027
        gamma_eff = gamma_0 - 1j * A_0 ** 2

        A = np.real(gamma_eff) - mu
        B = (nu + alpha / (sigma ** 2) + 2 * np.absolute(A_0) ** 2) + np.imag(gamma_eff)
        C = -(nu + alpha / (sigma ** 2) + 2 * np.absolute(A_0) ** 2) + np.imag(gamma_eff)
        D = -np.real(gamma_eff) - mu
        root = A ** 2 - 2 * A * D + 4 * B * C + D ** 2
        print(root[0])
        root = root.astype('complex128')
        omega_1 = 0.5 * (A + D + np.sqrt(root))
        omega_2 = 0.5 * (A + D - np.sqrt(root))
        OMEGAS_1.append(omega_1)
        OMEGAS_2.append(omega_2)

    ### PLOTS ###
    for i in range(len(OMEGAS_1)):
        plt.plot(phase, np.real(OMEGAS_1[i]), c=((gammas[i] - mu) / (gamma_max - mu), 0, 0), label=r"$\gamma_0=%.2f$" % (gammas[i],))
    plt.legend(loc=1)
    plt.xlim([-np.pi, np.pi])
    plt.grid(alpha=0.5)
    plt.xlabel("$\phi$", fontsize=20)
    plt.ylabel("$\\textrm{Re}(\lambda_+)$", fontsize=20)
    plt.tight_layout()
    plt.savefig('freq_real_+.png', dpi=300)
    plt.close()

    for i in range(len(OMEGAS_1)):
        plt.plot(phase, np.real(OMEGAS_2[i]), c=((gammas[i] - mu) / (gamma_max - mu), 0, 0), label=r"$\gamma_0=%.2f$" % (gammas[i],))
    plt.legend(loc=1)
    plt.xlim([-np.pi, np.pi])
    plt.grid(alpha=0.5)
    plt.xlabel("$\phi$", fontsize=20)
    plt.ylabel("$\\textrm{Re}(\lambda_-)$", fontsize=20)
    plt.tight_layout()
    plt.savefig('freq_real_-.png', dpi=300)
    plt.close()

    for i in range(len(OMEGAS_1)):
        plt.plot(phase, np.imag(OMEGAS_1[i]), c=((gammas[i] - mu) / (gamma_max - mu), 0, 0), label=r"$\gamma_0=%.2f$" % (gammas[i],))
    plt.legend(loc=1)
    plt.xlim([-np.pi, np.pi])
    plt.grid(alpha=0.5)
    plt.xlabel("$\phi$", fontsize=20)
    plt.ylabel("$\\textrm{Im}(\lambda_+)$", fontsize=20)
    plt.tight_layout()
    plt.savefig('freq_img_+.png', dpi=300)
    plt.close()

    for i in range(len(OMEGAS_1)):
        plt.plot(phase, np.imag(OMEGAS_2[i]), c=((gammas[i] - mu) / (gamma_max - mu), 0, 0), label=r"$\gamma_0=%.2f$" % (gammas[i],))
    plt.legend(loc=1)
    plt.xlim([-np.pi, np.pi])
    plt.grid(alpha=0.5)
    plt.xlabel("$\phi$", fontsize=20)
    plt.ylabel("$\\textrm{Im}(\lambda_-)$", fontsize=20)
    plt.tight_layout()
    plt.savefig('freq_img_-.png', dpi=300)
    plt.close()