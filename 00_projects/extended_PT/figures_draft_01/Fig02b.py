import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    directory = "C:/mnustes_science/simulation_data/FD/PDNLS_extended_PT/extras/dimensional/analysis"

    data = np.loadtxt(directory + '/data.txt', delimiter=',')
    data_dist = np.loadtxt(directory + '/data_dist.txt', delimiter=',')

    init = -9
    fin = -1

    ### FREQUENCY VS DISTANCE ###
    dist_full = data_dist[:, 0]
    freq_dist_full = 2 * np.pi * 1000 / data_dist[:, 1]
    freq_dist_err_full = 2 * np.pi * 1000 * np.abs(data_dist[:, 2]/data_dist[:, 1]**2)

    init = -20
    dist_bif = dist_full[init:]
    freq_dist_bif = freq_dist_full[init:]
    freq_dist_err_bif = freq_dist_err_full[init:]

    dist_01 = dist_full[1:21]
    freq_dist_01 = freq_dist_full[1:21]
    freq_dist_err_01 = freq_dist_err_full[1:21]

    dist_02 = dist_full[21:38]
    freq_dist_02 = freq_dist_full[21:38]
    freq_dist_err_02 = freq_dist_err_full[21:38]

    def F(x, A, c, b):
        return A * (((np.abs((x - c)) ) ** b))
    def G(x, a, b):
        return a * x + b
    def H(x, A, c):
        return A * (((np.abs((x - c))) ** 0.5))

    popt, pcov = curve_fit(F, dist_bif, freq_dist_bif, bounds=[(0, 30.28, 0.4), (100, 31, 0.51)])
    popt01, pcov01 = curve_fit(G, dist_01, freq_dist_01)
    popt02, pcov02 = curve_fit(G, dist_02, freq_dist_02)
    dist_dense = np.arange(dist_full[init - 1], popt[1], 0.00001)

    #[A, dc, alpha] = popt
    [A, dc, alpha] = popt
    [a, b] = popt01
    [c, d] = popt02

    delta_d = dist_bif[1] - dist_bif[0]

    fig, ax1 = plt.subplots()
    ax1.set_xlabel('$d\ \\textrm{(mm)}$', fontsize=30)
    ax1.set_ylabel('$\Omega$', fontsize=30)
    ax1.tick_params(labelsize=20)
    ax1.scatter(dist_full, freq_dist_full, c="k", zorder=10)
    zeros = np.arange(30.3, 33, 0.25)
    ax1.scatter(zeros, 0 * zeros, c="k")
    #ax1.plot(dist_01, G(dist_01, *popt01) - 0.5 * 2 * np.pi, linewidth=2)
    ax1.plot(dist_02, G(dist_02, *popt02) - 0.5 * 2 * np.pi, linewidth=2, c="r", linestyle="--")
    ax1.hlines(0, dist_full[0] - delta_d, 33, colors="k")
    #ax1.set_xlim([dist_full[0] - delta_d, dist_full[-1] + 10 * delta_d])
    ax1.set_xlim([25, 32])
    ax1.set_ylim([-1, 60])
    ax1.grid(alpha=0.2, zorder=0)


    left, bottom, width, height = [0.55, 0.56, 0.36, 0.34]
    ax2 = fig.add_axes([left, bottom, width, height])
    ax2.plot(dist_dense, F(dist_dense, popt[0], popt[1], popt[2]), color="r")
    ax2.errorbar(dist_bif, freq_dist_bif, freq_dist_err_bif, marker='o', ls='', ecolor="k", mec='black', color="k")
    ax2.grid(alpha=0.2)
    ax2.tick_params(labelsize=14)
    ax2.hlines(0, dist_bif[0] - delta_d, dc + delta_d, colors="k")
    ax2.set_xlim([dist_bif[0] - delta_d, dc + delta_d])
    text_01 = '\n'.join((r'$\Omega_0=%.3f$' % (A,), r'$d_c=%.3f$' % (dc,), r'$\alpha=%.4f$' % (alpha,)))
    text_02 = "$\Omega = \Omega_0 (d - d_c)^{\\alpha}$"
    text_03 = "$\Omega =" + r"%.2f" % (a,) + "d + " + r"%.2f" % (b,) + "$"
    text_04 = "$\Omega =" + r"%.2f" % (c,) + "d + " + r"%.2f" % (d,) + "$"
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax2.text(0.52, 0.7, text_01, transform=ax1.transAxes, fontsize=14, verticalalignment='top', bbox=props)
    ax2.text(0.68, 0.9, text_02, transform=ax1.transAxes, fontsize=14, verticalalignment='top')
    #ax2.text(0.05, 0.85, text_03, transform=ax1.transAxes, fontsize=16, verticalalignment='top', rotation=-32)
    ax2.text(0.1, 0.6, text_04, transform=ax1.transAxes, fontsize=20, verticalalignment='top', rotation=-30)
    ax2.text(0.5, 0.98, "$\\times 10^{-3}$", transform=ax1.transAxes, fontsize=12, verticalalignment='top')
    plt.tight_layout()
    #plt.show()
    plt.savefig('Fig02b.png', dpi=300)
    plt.close()
