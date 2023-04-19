import matplotlib.pyplot as plt
import numpy as np

from functions import *
from back_process import *
from time_integrators import *
from directories_lyap import *

if __name__ == '__main__':
    disco = 'E:/'
    initial_dir_data = str(disco) + "mnustes_science/simulation_data/FD/PDNLS_chaosa" + subdir_branching
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    #H_points_max = np.loadtxt(directory + '/H_max.txt', delimiter=',')[1:]
    H_points = np.loadtxt(directory + '/H_ext.txt', delimiter=',')
    T_points = np.loadtxt(directory + '/T_ext.txt', delimiter=',')
    #T_points_min = np.loadtxt(directory + '/T_min.txt', delimiter=',')[1:]
    H_mean = np.loadtxt(directory + '/H_mean.txt', delimiter=',')
    T = np.loadtxt(directory + '/T.txt', delimiter=',')
    print(len(H_mean))
    print(len(T))
    N_bins = 15
    cuts = [-15, -16, -17, -18, -19, -20]
    colors = []
    n_color = 7
    cut_str = []
    for cut in cuts:
        cut_str.append("$\Delta H = " + str(cut) + "$")
    for i in range(n_color):
        color_i = np.random.rand(1, 3)
        colors.append(color_i)
    histograms = []
    for k in range(len(cuts)):
        print(cuts[k])
        delta_H = []
        T_int = []
        for i in range(1, len(H_points)):
            delta_H_i = H_points[i] - H_points[i - 1]
            if delta_H_i < cuts[k]:
                T_int.append(T_points[i - 1])
            delta_H.append(delta_H_i)
        #plt.vlines(T_int, 0, 30, zorder=0, colors=colors[k])
        delta_T = []
        np.savetxt(directory + '/T_intervalos.txt', T_int, delimiter=',')
        for i in range(1, len(T_int)):
            delta_T.append(T_int[i] - T_int[i - 1])
        delta_T = np.array(delta_T)
        histograms.append(delta_T)
        hist, bins = np.histogram(delta_T, N_bins)
        dts = []
        for i in range(1, len(hist) + 1):
            dts.append((bins[i] + bins[i - 1]) / 2)
        plt.scatter(dts, np.log(hist), label=cut_str[k])
        plt.plot(dts, np.log(hist))
    plt.legend()
    plt.xlabel('$\Delta t$', size='20')
    plt.xticks(fontsize=15)
    # plt.xlim([T[0], T[-1]])
    plt.ylabel('$\\textrm{ln}\left(\\textrm{PDF}\\right)$', size='20')
    plt.yticks(fontsize=15)
    plt.ylim([0, 6])
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig(directory + '/dt_histogram_log_' + str(N_bins) + 'bins.png', dpi=300)
    plt.close()


    plt.hist(histograms, N_bins, histtype='bar', label=cut_str)
    plt.legend()
    plt.xlabel('$\Delta t$', size='20')
    plt.xticks(fontsize=15)
    #plt.xlim([T[0], T[-1]])
    plt.ylabel('$\\textrm{PDF}(\Delta t)$', size='20')
    plt.yticks(fontsize=15)
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig(directory + '/dt_histogram_' + str(N_bins) + 'bins.png', dpi=300)
    plt.close()


