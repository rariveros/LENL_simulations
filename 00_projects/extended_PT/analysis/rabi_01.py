from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    disco = 'D:/'
    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FDa'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    print("#############   " + directory + "   #############")
    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    X = np.loadtxt(directory + '/X.txt', delimiter=',')
    T = np.loadtxt(directory + '/T.txt', delimiter=',')
    params = np.loadtxt(directory + '/parameters.txt', delimiter=',')

    Nx = len(X)
    Nt = len(T)

    ti, tf = int(0.5 * Nt), Nt
    Z_r = Z_r[ti:tf, :]
    Z_i = Z_i[ti:tf, :]
    T = T[ti:tf] - T[ti]
    dt = T[1] - T[0]

    Nx = len(X)
    Nt = len(T)

    Z = Z_r + 1j * Z_i
    Z_conj = Z_r - 1j * Z_i
    Z_modulo = np.absolute(Z)

    Z_mod_L = Z_modulo[:, 0:int(Nx / 2)]
    Z_mod_R = Z_modulo[:, int(Nx / 2) + 1:]
    X_L = X[0:int(Nx / 2)]
    X_R = X[int(Nx / 2) + 1:]

    N_L = integrate.simpson(Z_mod_L ** 2, X_L)
    N_R = integrate.simpson(Z_mod_R ** 2, X_R)
    N = integrate.simpson(Z_modulo ** 2, X)

    plt.plot(T, N_L)
    plt.plot(T, N_R)
    plt.show()
    plt.close()