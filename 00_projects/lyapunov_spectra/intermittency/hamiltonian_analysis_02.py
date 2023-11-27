import matplotlib.pyplot as plt
import numpy as np

from functions import *
from back_process import *
from time_integrators import *
from directories_lyap import *
from scipy.interpolate import splrep, BSpline
from scipy.interpolate import splev, splrep
from scipy import interpolate
import warnings
warnings.simplefilter('ignore', np.RankWarning)

if __name__ == '__main__':
    disco = 'E:/'
    initial_dir_data = str(disco) + "mnustes_science/simulation_data/FD/PDNLS_chaosa" + subdir_branching
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    H_points = np.loadtxt(directory + '/H_ext.txt', delimiter=',')
    T_points = np.loadtxt(directory + '/T_ext.txt', delimiter=',')
    H_mean = np.loadtxt(directory + '/H_mean.txt', delimiter=',')
    print("Nt: " + str(len(H_mean)))
    N_bins = 60#0.5 * np.amax(H_mean)
    cuts = 1 * [0.5 * np.amax(H_mean)]#np.arange(30, 40, 2)
    colors = []
    n_color = 7
    cut_str = []
    HIST = []
    histograms = []
    for k in range(len(cuts)):
        delta_H = []
        T_int = []
        for i in range(1, len(H_points)):
            delta_H_i = H_points[i] - H_points[i - 1]
            if delta_H_i > cuts[k]:
                T_int.append((T_points[i - 1] + T_points[i]) / 2)
            delta_H.append(delta_H_i)
        delta_T = []
        for i in range(1, len(T_int)):
            delta_T.append(T_int[i] - T_int[i - 1])
        delta_T = np.array(delta_T)
        histograms.append(delta_T)
        hist, bins = np.histogram(delta_T, N_bins)
        dts = []
        for i in range(1, len(hist) + 1):
            dts.append((bins[i] + bins[i - 1]) / 2)
        hist = hist.astype('float')
        plt.scatter(dts, np.log(hist))
        HIST.append(hist)
    HIST = np.array(HIST)
    print(delta_T)
    print(bins)

    plt.legend()
    plt.xlabel('$\Delta t$', size='20')
    plt.xticks(fontsize=15)
    #plt.xlim([90, 600])
    plt.ylabel('$\\textrm{ln}\left(\\textrm{PDF}\\right)$', size='20')
    plt.yticks(fontsize=15)
    #plt.ylim([0, 4])
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


