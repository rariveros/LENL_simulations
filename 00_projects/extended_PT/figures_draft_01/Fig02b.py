import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    directory = "C:/mnustes_science/simulation_data/FD/PDNLS_extended_PT/extras/dimensional/analysis"

    #data = np.loadtxt(directory + '/data.txt', delimiter=',')
    data_dist = np.loadtxt(directory + '/data_dist.txt', delimiter=',')

    init = 0
    fin = -1

    ### FREQUENCY VS DISTANCE ###
    dist_full = data_dist[:fin, 0]
    freq_dist_full = 2 * np.pi * 1000 / data_dist[:fin, 1]
    print(dist_full)
    print(freq_dist_full)
    freq_dist_err_full = 2 * np.pi * 1000 * np.abs(data_dist[:fin, 2]/data_dist[:fin, 1]**2)

    init = -20
    dist_bif = dist_full[init:]
    freq_dist_bif = freq_dist_full[init:]
    freq_dist_err_bif = freq_dist_err_full[init:]

    dist_01 = dist_full[1:21]
    freq_dist_01 = freq_dist_full[1:21]
    freq_dist_err_01 = freq_dist_err_full[1:21]

    dist_02 = dist_full[0:46]
    freq_dist_02 = freq_dist_full[0:46]
    freq_dist_err_02 = freq_dist_err_full[0:46]

    #init = -10
    #dist_bif = dist_full[init:]
    #freq_dist_bif = freq_dist_full[init:]
    #freq_dist_err_bif = freq_dist_err_full[init:]

    #dist_01 = dist_full[:11]
    #freq_dist_01 = freq_dist_full[:11]
    #freq_dist_err_01 = freq_dist_err_full[:11]

    #dist_02 = dist_full[0:23]
    #freq_dist_02 = freq_dist_full[0:23]
    #freq_dist_err_02 = freq_dist_err_full[0:23]

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
    ax1.set_xlabel('$d\ \\textrm{(mm)}$', fontsize=25)
    ax1.set_ylabel('$\Omega \\times 10^{-3}$', fontsize=25)
    ax1.tick_params(labelsize=20)
    ax1.plot(np.arange(0, 33, 1), G(np.arange(0, 33, 1), *popt02), linewidth=2, c="r", zorder=0)
    ax1.scatter(dist_full, freq_dist_full, c="k", s=20, zorder=10)
    #props = dict(boxstyle='round,pad=0.1', facecolor='white', alpha=0.75)
    ax1.text(16.5, 125.5, "$\\textbf{(e)}$", fontsize=23, rotation=0, zorder=13)
    zeros_01 = np.arange(30.3, 35, 0.25)
    zeros_02 = np.arange(14, 18.5, 0.25)
    ax1.scatter(zeros_01, 0 * zeros_01, c="k", s=20, zorder=10)
    ax1.scatter(zeros_02, 0 * zeros_02, c="k", s=20, zorder=10)
    #ax1.plot(dist_01, G(dist_01, *popt01) - 0.5 * 2 * np.pi, linewidth=2)
    ax1.hlines(0, 0, 33, colors="k")
    ax1.vlines(18.41, 0, 150, colors="k", linestyles="--")
    ax1.vlines(30.2930, 0, 150, colors="k", linestyles="--")
    #ax1.set_xlim([dist_full[0] - delta_d, dist_full[-1] + 10 * delta_d])
    ax1.set_xlim([16, 32.5])
    ax1.set_ylim([-1, 140])
    ax1.grid(alpha=0.2, zorder=0)


    left, bottom, width, height = [0.6, 0.55, 0.34, 0.34]
    ax2 = fig.add_axes([left, bottom, width, height])
    ax2.plot(dist_dense, F(dist_dense, popt[0], popt[1], popt[2]), color="r")
    ax2.errorbar(dist_bif, freq_dist_bif, freq_dist_err_bif, marker='o', ls='', ecolor="k", mec='black', color="k")
    ax2.grid(alpha=0.2)
    ax2.tick_params(labelsize=14)
    ax2.hlines(0, dist_bif[0] - delta_d, dc + delta_d, colors="k")
    ax2.set_xlim([dist_bif[0] - delta_d, dc + delta_d])
    ax2.set_xlabel('$d\ \\textrm{(mm)}$', fontsize=17)
    ax2.set_ylabel('$\Omega \\times 10^{-3}$', fontsize=17)
    text_01 = '\n'.join((r'$\Omega_0=%.3f$' % (A,), r'$d_c=%.3f$' % (dc,), r'$\alpha=%.4f$' % (alpha,)))
    text_02 = "$\Omega = \Omega_0 (d - d_c)^{\\alpha}$"
    text_03 = "$\Omega =" + r"%.2f" % (a,) + "d + " + r"%.2f" % (b,) + "$"
    text_04 = "$\Omega(d) =" + r"%.2f" % (c,) + "d + " + r"%.2f" % (d,) + "$"
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax2.text(0.58, 0.7, text_01, transform=ax1.transAxes, fontsize=14, verticalalignment='top', bbox=props)
    ax2.text(0.69, 0.87, text_02, transform=ax1.transAxes, fontsize=14, verticalalignment='top')
    #ax2.text(0.05, 0.85, text_03, transform=ax1.transAxes, fontsize=16, verticalalignment='top', rotation=-32)
    ax2.text(0.29, 0.53, text_04, transform=ax1.transAxes, fontsize=20, verticalalignment='top', rotation=-33, color="r")
    #ax2.text(0.5, 0.98, "$\\times 10^{-3}$", transform=ax1.transAxes, fontsize=12, verticalalignment='top')
    ax2.text(0.10, -0.05, "$d_{RO}$", transform=ax1.transAxes, fontsize=25, verticalalignment='top')
    ax2.text(0.75, -0.05, "$d_c$", transform=ax1.transAxes, fontsize=25, verticalalignment='top')
    plt.tight_layout()
    #plt.show()
    plt.savefig('Fig02b.png', dpi=300)
    plt.close()
