import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':

    disco = 'C:/'
    initial_dir_data = str(disco) + ':/mnustes_science/simulation_data'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    X = np.loadtxt(directory + '/X.txt', delimiter=',')
    T = np.loadtxt(directory + '/T.txt', delimiter=',')

    Nx = len(X)

    U_complex = Z_r + 1j * Z_i
    X = X[int(0.1 * Nx): int(0.9 * Nx)]
    # Definiendo variables finales
    modulo_light = np.absolute(U_complex[-1, int(0.1 * Nx):int(0.9 * Nx)])
    arg = np.angle(U_complex[-1, int(0.1 * Nx):int(0.9 * Nx)])
    plt.plot(X, np.unwrap(arg, period=2 * np.pi) - 5.4, zorder=1, label="$\\theta(x)$", color="r")
    plt.plot(X, modulo_light, label="$R(x)$", color="b")
    plt.plot(X, np.exp(-X**2/(2*6**2)), label="$\gamma(x)$", color="k")
    plt.xlabel("$x$", fontsize=20)
    #plt.ylabel("$\\textrm{arg}(A)$", fontsize=20)
    plt.legend(fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.grid(alpha=0.3, zorder=0)
    plt.show()