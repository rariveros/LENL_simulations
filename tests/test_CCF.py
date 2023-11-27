import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    x_grid = np.arange(-20, 20, 0.01)
    Nx = len(x_grid)
    k = 2.3
    sigma = 2
    phi = np.cos(k * x_grid) * np.exp(- (x_grid - 2) ** 2 / (2 * sigma ** 2)) + np.cos(2 * k * x_grid) * np.exp(- (x_grid + 2) ** 2 / (2 * sigma ** 2))
    def CCF_wavelength_cos(k_min, k_max, dk):
        k_test = np.arange(k_min, k_max, dk)
        CCFs = []
        CCF_max = []
        for i in range(len(k_test)):
            phi_test_i = np.cos(k_test[i] * x_grid)
            CCF_i = np.correlate(phi, phi_test_i, "same")
            index = np.argmax(CCF_i)
            CCF_max_i = CCF_i[index]
            CCFs.append(CCF_i)
            CCF_max.append(CCF_max_i)
        CCFs = np.array(CCFs)
        CCF_max = np.array(CCF_max)
        ind = np.unravel_index(np.argmax(CCFs, axis=None), CCFs.shape)

        #graficar máximos de cada CCF en un grafico CCF_max v/s k
        fig, ax = plt.subplots()

        textstr = "$k_{\\textrm{max}}=" + str('%.2f' % k_test[ind[0]]) + "$"
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.75, 0.95, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)

        plt.plot(k_test, CCF_max, color="k")
        plt.scatter(k_test[ind[0]], CCF_max[ind[0]], c="r", zorder=10)

        plt.xlim(k_test[0], k_test[-1])
        plt.ylim([0, CCF_max[ind[0]] * 1.1])
        plt.grid(alpha=0.5)
        plt.tight_layout()
        plt.savefig("test_CCF.png", dpi=300)
        return k_test, CCF_max, k_test[ind[0]]


    k_test, CCF_max_k, k_max = CCF_wavelength_cos(0, 10, 0.001)
