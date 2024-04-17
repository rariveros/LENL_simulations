import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    disco = 'D:/'
    initial_dir_data = str(disco) + 'Users/mnustes_science/PT_fluids/mnustes_science/simulation_data'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    save_directory = directory + "/analysis"


    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    X = np.loadtxt(directory + '/X.txt', delimiter=',')
    T = np.loadtxt(directory + '/T.txt', delimiter=',')
    params = np.loadtxt(directory + '/parameters.txt', delimiter=',')
    gamma_real = np.loadtxt(directory + '/gamma_real.txt', delimiter=',')

    Nx = len(X)
    Nt = len(T)

    J_L = np.argmin(np.abs(X + params[3] / 2))
    J_R = np.argmin(np.abs(X - params[3] / 2))

    ti, tf = int(0.2 * Nt), Nt
    Z_r = Z_r[ti:tf, :]
    Z_i = Z_i[ti:tf, :]
    T = T[ti:tf] - T[ti]
    dt = T[1] - T[0]
    Nx = len(X)
    Nt = len(T)

    plt.plot(T, Z_r[:, J_R])
    plt.plot(T, Z_i[:, J_L])
    plt.show()
    plt.close()

    Z = Z_r + 1j * Z_i
    Z_last = Z[-1, :]
    N = 1
    theta_i = 0
    theta_f = 2 * np.pi
    thetas = [0]#np.arange(theta_i, theta_f, (theta_f - theta_i) / N)
    argmaxs = argrelextrema(np.abs(Z_last), np.greater, axis=0)
    envelope_R = np.abs(hilbert(np.real(Z_last)))
    envelope_I = np.abs(hilbert(np.imag(Z_last)))
    plt.plot(X, np.real(Z_last), color="b", label="$\psi_R$")
    plt.plot(X, np.imag(Z_last), color="r", label="$\psi_I$")
    plt.plot(X, gamma_real, color="k", label="$\gamma(x)$", linestyle="--")
    #plt.plot(X, envelope_I, color="r", label="$C_I(x)$", linestyle="--")
    #plt.plot(X, np.abs(Z_last), color="k", label="$|\psi|$")
    #plt.plot(X, np.sqrt(envelope_R ** 2 + envelope_I ** 2), color="k", label="$|C|)$", linestyle="--")

    def bigaussian(x, A1, a1, x01, A2, a2, x02):
        return A1 * np.exp(-(x - x01) ** 2 /(2 * a1 ** 2)) + A2 * np.exp(-(x - x02) ** 2 / (2 * a2 ** 2))


    #popt, pcov = curve_fit(bigaussian, X, np.sqrt(envelope_R ** 2 + envelope_I ** 2), bounds=([0.2, 10, -np.inf, 0.2, 10, 0], [np.inf, 40, 0, np.inf, 40, np.inf]))
    #plt.plot(X, bigaussian(X, *popt), color="k")
    #print(*popt)

    #if not os.path.exists(save_directory):
    #    os.makedirs(save_directory)

    #for i in range(len(argmaxs)):
    #    plt.scatter(X[argmaxs[i]], np.abs(Z_last[argmaxs[i]]))

    #for theta in thetas:
    #    Z_transf = Z_last
    #    plt.plot(X, np.real(Z_transf), color="b", label="$\psi_R$")
    #    plt.plot(X, np.imag(Z_transf), color="r", label="$\psi_I$")
    #    plt.plot(X, np.abs(Z_transf), color="k", label="$|\psi|$")

    plt.legend()
    plt.show()
    plt.close()
