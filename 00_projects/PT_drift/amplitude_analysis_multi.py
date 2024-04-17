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
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    directories = [name for name in os.listdir(working_directory) if os.path.isdir(os.path.join(working_directory, name))]
    for directory in directories:
        working_directory_i = working_directory + "/" + directory + "/sigma=16.000/gamma=0.2400"
        directories_i = [name for name in os.listdir(working_directory_i) if os.path.isdir(os.path.join(working_directory_i, name))]
        print(directories_i)
        values_01 = []
        values_02 = []

        for directory_j in directories_i:
            Z_r = np.loadtxt(working_directory_i + "/" + directory_j + '/field_real.txt', delimiter=',')
            Z_i = np.loadtxt(working_directory_i + "/" + directory_j + '/field_img.txt', delimiter=',')
            X = np.loadtxt(working_directory_i + "/" + directory_j + '/X.txt', delimiter=',')
            T = np.loadtxt(working_directory_i + "/" + directory_j + '/T.txt', delimiter=',')
            params = np.loadtxt(working_directory_i + "/" + directory_j + '/parameters.txt', delimiter=',')

            Nx = len(X)
            Nt = len(T)

            ti, tf = int(0.2 * Nt), Nt
            Z_r = Z_r[ti:tf, :]
            Z_i = Z_i[ti:tf, :]
            T = T[ti:tf] - T[ti]
            dt = T[1] - T[0]
            Nx = len(X)
            Nt = len(T)

            Z = Z_r + 1j * Z_i
            envelope_abs = []
            for i in range(Nt):
                envelope_R = np.abs(hilbert(np.real(Z[i, :])))
                envelope_I = np.abs(hilbert(np.imag(Z[i, :])))
                envelope_abs_i = np.sqrt(envelope_R ** 2 + envelope_I ** 2)
                envelope_abs.append(envelope_abs_i)
            envelope_abs_mean = np.mean(np.array(envelope_abs), axis=0)

            def bigaussian(x, A1, a1, x01, A2, a2, x02):
                return A1 * np.exp(-(x - x01) ** 2 /(2 * a1 ** 2)) + A2 * np.exp(-(x - x02) ** 2 / (2 * a2 ** 2))

            popt, pcov = curve_fit(bigaussian, X, envelope_abs_mean, bounds=([0.2, 10, -np.inf, 0.2, 10, 0], [np.inf, 40, 0, np.inf, 40, np.inf]))
            A1 = popt[0]
            A2 = popt[3]
            if A1 > A2:
                A1 = popt[0]
                A2 = popt[3]
            else:
                A2 = popt[0]
                A1 = popt[3]
            values_01.append([A1, params[3]])
            values_02.append([A2, params[3]])
        values_01 = np.array(values_01)
        values_02 = np.array(values_02)
        plt.scatter(values_01[:, 1], (values_01[:, 0] - values_02[:, 0]) / (values_01[:, 0] + values_02[:, 0]))
        #plt.scatter(values_02[:, 1], values_02[:, 0])
    plt.show()


