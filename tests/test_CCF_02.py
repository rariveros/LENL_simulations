import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    disco = 'D:/'
    initial_dir_data = str(disco) + 'Users/mnustes_science/PT_fluids/mnustes_science/simulation_data'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    X = np.loadtxt(directory + '/X.txt', delimiter=',')
    T = np.loadtxt(directory + '/T.txt', delimiter=',')
    Nx = len(X)
    Nt = len(T)
    ti = int(0.5 * Nt)
    tf = -1

    Z_r = Z_r[ti:tf, :]
    Z_i = Z_i[ti:tf, :]
    T = T[ti:tf]
    Nt = len(T)
    Z = Z_r + 1j * Z_i
    modulo = np.absolute(Z)

    def CCF_wavelength_cos(x_grid, phi, k_min, k_max, dk, plot):
        k_test = np.arange(k_min, k_max, dk)
        CCFs = []
        CCF_max = []
        for i in range(len(k_test)):
            phi_test_i = np.absolute(np.cos((2 * np.pi / k_test[i]) * x_grid))
            CCF_i = np.correlate(phi, phi_test_i, "same")
            index = np.argmax(CCF_i)
            CCF_max_i = CCF_i[index]
            CCFs.append(CCF_i)
            CCF_max.append(CCF_max_i)
        CCFs = np.array(CCFs)
        CCF_max = np.array(CCF_max)
        ind = np.unravel_index(np.argmax(CCFs, axis=None), CCFs.shape)
        if plot == "SI":
            fig, ax = plt.subplots()

            textstr = "$\lambda_{\\textrm{max}}=" + str('%.2f' % k_test[ind[0]]) + "$"
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            ax.text(0.75, 0.95, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)

            plt.plot(k_test, CCF_max, color="k")
            plt.scatter(k_test[ind[0]], CCF_max[ind[0]], c="r", zorder=10)

            plt.xlim(k_test[0], k_test[-1])
            plt.ylim([0, CCF_max[ind[0]] * 1.1])
            plt.grid(alpha=0.5)
            plt.tight_layout()
            plt.savefig("test_CCF.png", dpi=300)
        return k_test, CCF_max, ind[0]


    x_min = int(3 * Nx / 8)
    x_center = int(4 * Nx / 8)
    x_max = int(5 * Nx / 8)
    k_test_01 = []
    CCF_max_k_01 = []
    ind_max_01 = []
    k_test_02 = []
    CCF_max_k_02 = []
    ind_max_02 = []

    for i in range(Nt):
        k_test_01_i, CCF_max_k_01_i, ind_max_01_i = CCF_wavelength_cos(X[x_min:x_max], Z_i[i, x_min:x_center], 10, 120, 0.01, plot="NO")
        k_test_02_i, CCF_max_k_02_i, ind_max_02_i = CCF_wavelength_cos(X[x_min:x_max], Z_r[i, x_center:x_max], 10, 120, 0.01, plot="NO")
        k_test_01.append(k_test_01_i)
        CCF_max_k_01.append(CCF_max_k_01_i)
        ind_max_01.append(ind_max_01_i)
        k_test_02.append(k_test_02_i)
        CCF_max_k_02.append(CCF_max_k_02_i)
        ind_max_02.append(ind_max_02_i)

    k_test_01 = np.array(k_test_01)
    CCF_max_k_01 = np.mean(np.array(CCF_max_k_01), axis=0)
    ind_max_01 = int(np.mean(np.array(ind_max_01)))
    k_test_02 = np.array(k_test_02)
    CCF_max_k_02 = np.mean(np.array(CCF_max_k_02), axis=0)
    ind_max_02 = int(np.mean(np.array(ind_max_02)))

    fig, ax = plt.subplots()

    #textstr = "$\lambda^{L}_{\\textrm{max}}=" + str('%.2f' % k_test_01[ind_max_01]) + "$\n" \
    #                                                                             + "$\lambda^{R}_{\\textrm{max}}=" + str('%.2f' % k_test_02[ind_max_02]) +"$"

    #props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    #ax.text(0.75, 0.95, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)

    plt.plot(k_test_01_i, CCF_max_k_01, color="k")
    plt.plot(k_test_02_i, CCF_max_k_02, color="k")
    plt.scatter(k_test_01_i[ind_max_01], CCF_max_k_01[ind_max_01], c="r", zorder=10)
    plt.scatter(k_test_02_i[ind_max_02], CCF_max_k_02[ind_max_02], c="r", zorder=10)

    plt.xlim(10, 120)
    plt.grid(alpha=0.5)
    plt.tight_layout()
    plt.savefig("test_CCF.png", dpi=300)
