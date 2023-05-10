from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    phase = np.arange(0, 2 * np.pi, 0.01)
    gammas = np.arange(0.1, 0.18, 0.01)
    for gamma_i in gammas:
        gamma_0 = gamma_i
        mu = 0.1
        amp_i = ((2 * mu / 9) * (gamma_i - mu)) ** 0.25
        A_0 = amp_i * np.exp(1j * phase)
        nu = 0.027
        gamma_eff = gamma_0 - 1j * A_0 ** 2

        A = np.real(gamma_eff) - mu
        B = (nu + 2 * np.abs(A_0) ** 2) + np.imag(gamma_eff)
        D = -(nu + 2 * np.abs(A_0) ** 2) + np.imag(gamma_eff)
        C = -np.real(gamma_eff) - mu

        omega_1 = 0.5 * (((A - C) - 1j * (B + D)) + np.sqrt((-(A - C) + 1j * (B + C)) ** 2 - 8 * 1j * (B*C-A*D)))
        omega_2 = 0.5 * (((A - C) - 1j * (B + D)) - np.sqrt((-(A - C) + 1j * (B + C)) ** 2 - 8 * 1j * (B*C-A*D)))
        plt.plot(phase, np.real(omega_1), c=((gamma_i - mu) / (0.18 - mu), 0, 0), label=r"$\gamma_0=%.2f$" % (gamma_i,))
    plt.legend(loc=1)
    plt.xlim([0, 2*np.pi])
    plt.grid(alpha=0.5)
    plt.xlabel("$\phi$", fontsize=20)
    plt.ylabel("$\\textrm{Re}(\Omega_+)$", fontsize=20)
    plt.savefig('freq_real_+.png', dpi=300)