from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    T = [13.6145, 14.5510, 16.0120, 17.9568]
    T_err = [0.1555, 0.4315, 0.6432, 1.3790]
    dist = [0.00, 0.40, 0.79, 1.19]
    plt.errorbar(dist, T, yerr=T_err, fmt='o')
    plt.show()
    plt.close()

    freq = 1 / np.array(T)
    freq_err = (1 / np.array(T) ** 2) * np.array(T_err)
    plt.errorbar(31.6 + np.array(dist),  1000 * freq, yerr=1000 * freq_err, fmt='o')
    plt.xlim(31, 38)
    plt.ylim(0, 130)
    plt.show()
    plt.close()