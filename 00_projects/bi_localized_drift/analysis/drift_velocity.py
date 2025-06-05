import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *
from scipy.signal import find_peaks
from scipy.stats import linregress
from scipy.signal import argrelmax

if __name__ == '__main__':
    frequencies = []
    disc = "aD:"
    initial_dir_data = str(disc) + 'mnustes_science/simulation_data/FD'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    X = np.loadtxt(directory + '/X.txt', delimiter=',')
    T = np.loadtxt(directory + '/T.txt', delimiter=',')
    ZR = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    ZI = np.loadtxt(directory+ '/field_img.txt', delimiter=',')
    params = np.loadtxt(directory+ '/parameters.txt', delimiter=',')
    Z = ZR + 1j * ZI
    Z_mod = np.abs(Z)

    sigma = params[4]
    d = params[3]
    JL = np.argmin(np.abs(X + d / 2))
    JR = np.argmin(np.abs(X - d / 2))
    J_sig = np.argmin(np.abs(X - sigma)) - np.argmin(np.abs(X))

    Nx = len(X)
    Nt = len(T)
    dx = X[1] - X[0]
    dt = T[1] - T[0]
    I0 = int(0.1 * Nt)

    Z_mod = Z_mod[I0:]
    T = T[I0:]
    Nt = len(T)

    T_measure = 8
    Xs = []
    Ts = []
    for j in range(T_measure):
        maximos_idx = argrelmax(Z_mod[:, JL + int(T_measure / 2) - j], order=5)[0]
        x_measure = [X[JL + int(T_measure / 2) - j]] * len(maximos_idx[1:-1])
        t_measure = T[maximos_idx[1:-1]]
        Xs.append(x_measure)
        Ts.append(t_measure)
    N_paths = len(maximos_idx[1:-1])
    print(len(np.array(Xs)[0, :]))
    print(N_paths)
    velocidades = []
    for i in range(N_paths):
        slope, intercept, r_value, p_value, std_err = linregress(np.array(Ts)[:, i], np.array(Xs)[:, i])
        velocidades.append(slope)
    v_mean = np.mean(velocidades)
    print(v_mean)
