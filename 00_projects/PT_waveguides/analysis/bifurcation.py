from functions import *
from back_process import *

if __name__ == '__main__':
    gamma = 10.0
    k = 1.0
    beta = 0.5

    # PARAMETROS INICIALES DE POTENCIA (N) Y PHASE INICIAL DE PRIMER OSCILADOR
    N = np.arange(0, 30, 0.02, dtype=complex)
    n = 0
    z_pp_0 = +((N + 2) ** 2 - (2 * gamma ** 2) / (k ** 2) + ((2 * gamma) / (k ** 2)) * (
                gamma ** 2 - 4 * k * (-1) ** n * (1 + N)) ** 0.5) ** 0.5
    z_mm_0 = -((N + 2) ** 2 - (2 * gamma ** 2) / (k ** 2) - ((2 * gamma) / (k ** 2)) * (
                gamma ** 2 - 4 * k * (-1) ** n * (1 + N)) ** 0.5) ** 0.5
    z_pm_0 = +((N + 2) ** 2 - (2 * gamma ** 2) / (k ** 2) - ((2 * gamma) / (k ** 2)) * (
                gamma ** 2 - 4 * k * (-1) ** n * (1 + N)) ** 0.5) ** 0.5
    z_mp_0 = -((N + 2) ** 2 - (2 * gamma ** 2) / (k ** 2) + ((2 * gamma) / (k ** 2)) * (
                gamma ** 2 - 4 * k * (-1) ** n * (1 + N)) ** 0.5) ** 0.5

    n = 1
    z_pp_1 = +((N + 2) ** 2 - (2 * gamma ** 2) / (k ** 2) + ((2 * gamma) / (k ** 2)) * (
                gamma ** 2 - 4 * k * (-1) ** n * (1 + N)) ** 0.5) ** 0.5
    z_mm_1 = -((N + 2) ** 2 - (2 * gamma ** 2) / (k ** 2) - ((2 * gamma) / (k ** 2)) * (
                gamma ** 2 - 4 * k * (-1) ** n * (1 + N)) ** 0.5) ** 0.5
    z_pm_1 = +((N + 2) ** 2 - (2 * gamma ** 2) / (k ** 2) - ((2 * gamma) / (k ** 2)) * (
                gamma ** 2 - 4 * k * (-1) ** n * (1 + N)) ** 0.5) ** 0.5
    z_mp_1 = -((N + 2) ** 2 - (2 * gamma ** 2) / (k ** 2) + ((2 * gamma) / (k ** 2)) * (
                gamma ** 2 - 4 * k * (-1) ** n * (1 + N)) ** 0.5) ** 0.5

    z_stable = np.zeros(len(N))
    z_unstable = np.zeros(len(N))

    tol = 1e-3
    cut_01, cut_02 = 0.25, 15.74

    z_pp_0 = np.where(np.abs(np.imag(z_pp_0)) > tol, np.nan, z_pp_0.real)
    z_mm_0 = np.where(np.abs(np.imag(z_mm_0)) > tol, np.nan, z_mm_0.real)
    z_pm_0 = np.where(np.abs(np.imag(z_pm_0)) > tol, np.nan, z_pm_0.real)
    z_mp_0 = np.where(np.abs(np.imag(z_mp_0)) > tol, np.nan, z_mp_0.real)

    z_pp_1 = np.where(np.abs(np.imag(z_pp_1)) > tol, np.nan, z_pp_1.real)
    z_mm_1 = np.where(np.abs(np.imag(z_mm_1)) > tol, np.nan, z_mm_1.real)
    z_pm_1 = np.where(np.abs(np.imag(z_pm_1)) > tol, np.nan, z_pm_1.real)
    z_mp_1 = np.where(np.abs(np.imag(z_mp_1)) > tol, np.nan, z_mp_1.real)

    # stable: nan inside [cut_01, cut_02]
    z_stable = np.where((N >= cut_01) & (N <= cut_02), np.nan, z_stable)

    # unstable: nan outside [cut_01, cut_02]
    z_unstable = np.where((N < cut_01) | (N > cut_02), np.nan, z_unstable)

    fig, ax = plt.subplots(1, 1, figsize=(8, 4))

    ax.plot(N, np.real(z_pp_0), color="k", lw=4)
    ax.plot(N, np.real(z_mm_0), color="k", lw=4, ls="--")
    ax.plot(N, np.real(z_pm_0), color="k", lw=4, ls="--")
    ax.plot(N, np.real(z_mp_0), color="k", lw=4)
    ax.plot(N, z_unstable, color="k", lw=4, ls="--")
    ax.plot(N, z_stable, color="k", lw=4)
    ax.vlines([14, 20, 26], -25, 25, colors="r", linewidths=2)

    ax.set_xlabel('$N$', fontsize=28)
    ax.set_ylabel('$z$', fontsize=28)
    ax.tick_params(axis="both", direction="in", labelsize=23)
    ax.set_xlim([0, 30])
    ax.set_ylim([-25, 25])

    plt.tight_layout()
    plt.savefig("bifurcation.png", dpi=300)