from functions import *
from back_process import *
from time_integrators import *

def H(x, A, c):
    return A * ((np.abs((x - c))) ** 0.5)

def F(x, A, c):
    return A * np.real(((c - x + 0j) ** 0.25))

if __name__ == '__main__':
    working_directory = "D:/mnustes_science/experimental_data/soliton_control/espacio_parametros/gamma_cte1"

    centers = np.loadtxt(working_directory + '/c_promedios_1.txt', delimiter=',')
    centers_std = np.loadtxt(working_directory + '/std_1.txt', delimiter=',')
    nus = np.loadtxt(working_directory + '/nus_1.txt', delimiter=',')

    fig, ax = plt.subplots(1, 1, figsize=(4.5, 2.6), dpi=300)
    ax.errorbar(nus[:-1], -centers[:-1], centers_std[:-1], marker='o', ls='', ecolor="k", mec='black', color="k", ms=8, zorder=10)
    ax.hlines(0, -0.11, -0.06, colors="k", lw=0.8)
    ax.set_xlim([-0.107, -0.064])
    ax.set_xlabel(r'$\nu$', size=17)
    ax.set_ylim([-2.2, 14.0])
    ax.set_ylabel(r'$x_s\ \textrm{(mm)}$', size=17)
    ax.tick_params(labelsize=15)

    ax.grid(linestyle='--', alpha=0.05, color='k')

    plt.tight_layout(rect=[0, 0.03, 1, 1])
    plt.savefig("Fig01c.png", dpi=300)
    plt.close()

