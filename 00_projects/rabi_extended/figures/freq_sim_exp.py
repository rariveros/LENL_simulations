from functions import *
from back_process import *
from time_integrators import *

def H(x, A, c):
    return A * ((np.abs((x - c))) ** 0.5)

def F(x, A, c):
    return A * np.real(((c - x + 0j) ** 0.5))

if __name__ == '__main__':
    working_directory = "D:/mnustes_science/experimental_data/rabi_04_16"
    sim_directory = working_directory + "/simulations"
    exp_directory = working_directory + "/experimental"

    dist_sim = np.loadtxt(sim_directory + '/dists.txt', delimiter=',')
    freq_sim = 1000 * np.loadtxt(sim_directory + '/frequencies.txt', delimiter=',')
    freq_std = 1000 * np.loadtxt(sim_directory + '/freq_std.txt', delimiter=',')

    dist_exp = np.loadtxt(exp_directory + '/dists.txt', delimiter=',')
    freqL_exp = 1000 * np.loadtxt(exp_directory + '/frequencies_L.txt', delimiter=',')
    freqR_exp = 1000 * np.loadtxt(exp_directory + '/frequencies_R.txt', delimiter=',')
    freqL_std_exp = 1000 * np.loadtxt(exp_directory + '/frequencies_L_STD.txt', delimiter=',')
    freqR_std_exp = 1000 * np.loadtxt(exp_directory + '/frequencies_R_STD.txt', delimiter=',')

    ni_fit = -21
    nf_fit = -10
    dist_dense = np.arange(42.9, 44.7, 0.004)

    #plt.scatter(dist_sim[ni_fit:nf_fit], freq_sim[ni_fit:nf_fit])
    #plt.show()
    #plt.close()

    popt, pcov = curve_fit(F, dist_sim[ni_fit:nf_fit], freq_sim[ni_fit:nf_fit], bounds=[(0, 43.5), (100, 44.5)])

    n1 = 4
    n2 = -25
    n3 = -10

    fig, ax1 = plt.subplots()
    ax1.set_xlabel('$d\ \\textrm{(mm)}$', fontsize=25)
    ax1.set_ylabel('$\Omega\ \\textrm{(mHz)}$', fontsize=25)
    ax1.tick_params(labelsize=20)
    ax1.plot(dist_sim[:n1], freq_sim[:n1], linewidth=3, c="r", zorder=5, alpha=0.5)
    ax1.scatter(dist_sim[:n1], freq_sim[:n1], s=70, marker="^", c="r", edgecolor="k", lw=0.8, zorder=5)

    ax1.plot(dist_sim[n1:], freq_sim[n1:], linewidth=3, c="r", zorder=5, alpha=0.5)
    ax1.scatter(dist_sim[n1:], freq_sim[n1:], s=70, marker="^", c="r", edgecolor="k", lw=0.8, zorder=5)

    ax1.errorbar(dist_exp[:-1], (freqR_exp[:-1] + freqL_exp[:-1]) / 2, (freqR_std_exp[:-1] + freqL_std_exp[:-1]) / 2, marker='o', ls='', ecolor="k", mec='black', color="k", ms=8, zorder=10)
    ax1.hlines(0, 30, 56, colors="k", zorder=1)
    ax1.vlines(31.75, 0, 150, colors="k", linestyles="--")
    ax1.vlines(popt[1], 0, 150, colors="k", linestyles="--", alpha=0.7)
    ax1.set_xlim([30, 51])
    ax1.set_ylim([-3, 130])
    ax1.grid(alpha=0.2, zorder=0)

    left, bottom, width, height = [0.625, 0.635, 0.325, 0.32]
    ax2 = fig.add_axes([left, bottom, width, height])
    ax2.plot(dist_dense, F(dist_dense, popt[0], popt[1]), linewidth=3, color="r", zorder=4, ls="--", alpha=0.5)
    ax2.scatter(dist_sim[n2:n3], freq_sim[n2:n3], s=70, marker="^", c="r", edgecolor="k", lw=0.8, zorder=5)
    ax2.grid(alpha=0.2)
    ax2.tick_params(labelsize=14)
    ax2.hlines(0, 40, 45, colors="k", zorder=1)
    delta_d = dist_sim[n2 + 2] - dist_sim[n2 + 1]
    ax2.set_xlim([dist_sim[n2] - delta_d, 44.3])
    ax2.set_ylim([-1, 30])
    ax2.set_xlabel('$d\ \\textrm{(mm)}$', fontsize=17)
    ax2.set_ylabel('$\Omega\ \\textrm{(mHz)}$', fontsize=17)
    A = popt[0]
    dc = popt[1]
    text_01 = '\n'.join((r'$\Omega_0=%.1f$' % (A,) + '\ \\textrm{mHz}', r'$d_c=%.1f$' % (dc,) + '\ \\textrm{mm}'))
    text_02 = "$\Omega = \Omega_0 (d_c - d)^{1/2}$"
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax2.text(0.59, 0.712, text_01, transform=ax1.transAxes, fontsize=13, verticalalignment='top', bbox=props)
    ax2.text(0.68, 0.96, text_02, transform=ax1.transAxes, fontsize=14, verticalalignment='top')
    plt.tight_layout()

    plt.tight_layout()
    #plt.show()
    plt.savefig(working_directory + '/sim_exp_comparation.png', dpi=300)
    plt.close()