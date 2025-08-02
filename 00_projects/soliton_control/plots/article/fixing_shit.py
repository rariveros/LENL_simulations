from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    working_directory = "D:/mnustes_science/experimental_data/soliton_control"
    nonshifted_soliton = working_directory + "/espacio_parametros/mov_400cu/a=11.29_f=13.80_1"

    Z_nonshifted = np.loadtxt(nonshifted_soliton + '/Z_strobo.txt', delimiter=',')
    X_nonshifted = np.loadtxt(nonshifted_soliton + '/X_mm.txt', delimiter=',')
    Nx = len(Z_nonshifted[0, :])
    print(int(Nx/2))
    X_pix = np.arange(0, Nx)
    plt.plot(X_pix, Z_nonshifted[0, :])

    plt.vlines(int(Nx/2) + 4, 0, np.max(Z_nonshifted[0, :]))
    plt.show()