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
    values_01 = []
    values_02 = []
    #[alpha, beta, gamma_0, dist, mu, nu, sigma]
    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    X = np.loadtxt(directory + '/X.txt', delimiter=',')
    T = np.loadtxt(directory + '/T.txt', delimiter=',')
    params = np.loadtxt(directory + '/parameters.txt', delimiter=',')

    Nx = len(X)
    Nt = len(T)

    ti, tf = int(0.2 * Nt), Nt
    tint = 10
    xi, xf = -100, 100
    ji, jf = np.argmin(np.abs(X - xi)), np.argmin(np.abs(X - xf))
    Z_r = Z_r[ti:tf:tint, ji:jf]
    Z_i = Z_i[ti:tf:tint, ji:jf]
    X = X[ji:jf]
    T = T[ti:tf:tint] - T[ti]
    dt = T[1] - T[0]
    Nx = len(X)
    Nt = len(T)

    sigma_i = params[6]
    Z = Z_r + 1j * Z_i
    Z_abs = np.abs(Z)

    envelope_abs = []
    argmaxs = argrelextrema(Z_abs[0, :], np.greater, axis=0)
    N_peaks = len(argmaxs[0])
    m = N_peaks
    for i in range(1, Nt):
        argmaxs = argrelextrema(Z_abs[i, :], np.greater, axis=0)
        N_peaks_new = len(argmaxs[0])
        if N_peaks_new > N_peaks:
            m = m + 1
        N_peaks = len(argmaxs[0])
        for arg in argmaxs[0]:
            plt.scatter(X[arg], T[i], c="k")
    print(m)
    M_max = m
    positions = np.empty((Nt, M_max))
    positions.fill(np.NaN)

    argmaxs = argrelextrema(Z_abs[0, :], np.greater, axis=0)
    N_peaks = len(argmaxs[0])
    k = N_peaks
    for i in range(1, Nt):
        argmaxs = argrelextrema(Z_abs[i, :], np.greater, axis=0)
        positions[i, -k] = argmaxs[0]
        N_peaks_new = len(argmaxs[0])
        if N_peaks_new > N_peaks:
            k = k + 1
        N_peaks = len(argmaxs[0])
        #for arg in argmaxs[0]:
            #plt.scatter(X[arg], T[i], c="k")
    print(m)
    plt.show()



