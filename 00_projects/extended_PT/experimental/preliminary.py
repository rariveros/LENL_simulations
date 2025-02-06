from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    T = [13.97, 13.25, 12.74, 12.97, 13.08, 12.60, 13.18, 12.84, 12.67, 12.51, 14.56, 13.16, 10000]
    T_err = [0.16, 0.13, 0.16, 0.11, 0.05, 0.15, 0.12, 0.15, 0.24, 0.22, 0.37, 0.57, 0]
    dist = [0.00, 0.40, 0.80, 1.20, 1.47, 1.89, 2.14, 2.43, 3.78, 4.89, 7.67, 10.27, 12.85]

    #plt.plot(31.6 + np.array(dist),  T, color="k", lw=2)
    plt.errorbar(
        31.6 + np.array(dist), T, yerr=T_err,
        fmt='o', color="k",
        ecolor="k",  # Color de las barras de error
        elinewidth=1.5,  # Grosor de la barra de error
        capsize=3,  # Tamaño de los caps
        capthick=1.5  # Grosor de los caps
    )
    plt.xlabel("$d\ \\textrm{(mm)}$", fontsize=20)
    plt.ylabel("$T\ \\textrm{(s)}$", fontsize=20)
    plt.xlim(31, 55)
    plt.ylim(-1, 20)
    plt.grid(alpha=0.2)
    plt.show()
    plt.close()

    freq = 1 / np.array(T)
    freq_err = (1 / np.array(T) ** 2) * np.array(T_err)
    #plt.plot(31.6 + np.array(dist), 1000 * freq, color="k", lw=2)
    plt.errorbar(31.6 + np.array(dist),  1000 * freq, yerr=1000 * freq_err,
        fmt='o', color="k",
        ecolor="k",  # Color de las barras de error
        elinewidth=1.5,  # Grosor de la barra de error
        capsize=3,  # Tamaño de los caps
        capthick=1.5  # Grosor de los caps
    )
    plt.xlabel("$d\ \\textrm{(mm)}$", fontsize=20)
    plt.ylabel("$f\ \\textrm{(mHz)}$", fontsize=20)
    plt.xlim(31, 55)
    plt.ylim(-4, 130)
    plt.grid(alpha=0.2)
    plt.show()
    plt.close()