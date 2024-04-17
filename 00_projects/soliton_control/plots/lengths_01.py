import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    disco = 'E:/'
    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FD'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    X = np.loadtxt(directory + '/X.txt', delimiter=',')
    T = np.loadtxt(directory + '/T.txt', delimiter=',')
    params = np.loadtxt(directory + '/parameters.txt', delimiter=',') #[alpha, beta, gamma_0, mu_0, nu]
    sigma = 6

    U_complex = Z_r + 1j * Z_i
    modulo = np.absolute(U_complex)
    arg = np.angle(U_complex)
    arg = (2 * np.pi + arg) * (arg < 0) + arg * (arg > 0)
    final_phase = arg[-1, :]
    final_amplitude = modulo[-1, :]
    forcing = params[2] * np.exp(- X ** 2 / (2 * sigma ** 2))
    delta_homo = - params[4] + np.sqrt(params[2] ** 2 + params[3] ** 2)
    def soliton(x, n, delta, x0):
        return n * np.sqrt(2 * delta) * (1 / (np.cosh(np.sqrt(delta) * (x - x0))))

    popt, pcov = curve_fit(soliton, X, final_amplitude, bounds=([0, 0, -30], [10, 10, 30]))
    print(np.sqrt(params[0]) / np.sqrt(popt[1]))
    print(2*np.sqrt(2 * np.log(2)) * sigma)
    plt.plot(X, soliton(X, *popt) / np.amax(final_amplitude), color="r", linestyle="--")
    plt.plot(X, final_phase / np.amax(final_phase))
    plt.plot(X, final_amplitude / np.amax(final_amplitude))
    plt.hlines(1 / 2, -50, 50)
    plt.plot(X, forcing / np.amax(forcing))
    plt.xlim([-50, 50])
    plt.show()