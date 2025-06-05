import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *
from scipy.signal import find_peaks
from scipy.stats import linregress
from scipy.signal import argrelmax


def mean_velocity(Z_mod, X, t_grid, T_measure, J_side):
    Xs = []
    Ts = []
    for j in range(T_measure):
        maximos_idx, _ = find_peaks(Z_mod[:, J_side + int(T_measure / 2) - j], prominence=0.0001)
        if j == 0:
            N_path = len(maximos_idx[1:-1])
            if N_path > 4:
                N_path = 4
            first_idx = maximos_idx[1]
        #print((maximos_idx[2] - maximos_idx[1]))
        #print(maximos_idx[1] - first_idx)
        if np.abs(maximos_idx[1] - first_idx) > 0.5 * (maximos_idx[2] - maximos_idx[1]):
            maxs_idx = maximos_idx[2:N_path + 2]
            #print("hola")
            #print(len(maxs_idx))
        else:
            maxs_idx = maximos_idx[1:N_path + 1]
            #print("chao")
            #print(len(maxs_idx))
        #plt.title(str(j))
        #plt.plot(T, Z_mod[:, J_side + int(T_measure / 2) - j])
        #plt.scatter(T[maxs_idx], Z_mod[maxs_idx, J_side + int(T_measure / 2) - j])
        #plt.show()
        x_measure = [X[J_side + int(T_measure / 2) - j]] * len(maxs_idx)
        t_measure = t_grid[maxs_idx]
        Xs.append(x_measure)
        Ts.append(t_measure)
    maximos_idx_0, _ = find_peaks(Z_mod[:, J_side + int(T_measure / 2)], prominence=0.0000001)
    amp_mean = np.mean(Z_mod[maximos_idx_0, J_side + int(T_measure / 2)])
    T_mean = t_grid[maximos_idx[2]] - t_grid[maximos_idx[1]]
    N_paths = len(maxs_idx)
    velocidades = []
    for i in range(N_paths):
        slope, intercept, r_value, p_value, std_err = linregress(np.array(Ts)[:, i], np.array(Xs)[:, i])
        velocidades.append(slope)
    v_mean = np.mean(velocidades)
    return v_mean, amp_mean, T_mean

if __name__ == '__main__':
    disc = "aD:"
    initial_dir_data = str(disc) + 'mnustes_science/simulation_data/FD'
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    directories = [name for name in os.listdir(working_directory) if os.path.isdir(os.path.join(working_directory, name))]
    velocities_R = []
    velocities_L = []
    amplitudes_R = []
    amplitudes_L = []
    periods_R = []
    periods_L = []
    gammas = []
    for directory_01 in directories:
        directory = working_directory + "/" + directory_01
        X = np.loadtxt(directory + '/dist=80.000/X.txt', delimiter=',')
        T = np.loadtxt(directory + '/dist=80.000/T.txt', delimiter=',')
        ZR = np.loadtxt(directory + '/dist=80.000/field_real.txt', delimiter=',', dtype=complex)
        ZI = np.loadtxt(directory+ '/dist=80.000/field_img.txt', delimiter=',', dtype=complex)
        params = np.loadtxt(directory+ '/dist=80.000/parameters.txt', delimiter=',')
        #params = np.loadtxt(directory + '/dist=80.000/parameters.txt', delimiter=',')
        Z = ZR + 1j * ZI
        Z_mod = np.abs(Z)

        sigma = params[4]
        d = params[3]
        gamma = params[2]
        print("### gamma = " + str(d) + " ###")

        JL = np.argmin(np.abs(X + d / 2))
        JR = np.argmin(np.abs(X - d / 2))
        J_sig = np.argmin(np.abs(X - sigma)) - np.argmin(np.abs(X))

        Nx = len(X)
        Nt = len(T)
        dx = X[1] - X[0]
        dt = T[1] - T[0]
        if np.abs(gamma - 1) < 0.001:
            I0 = int(0.1 * Nt)  # int(0.165 * Nt)
            print("hohola")
        else:
            I0 = int(0.1 * Nt)#int(0.165 * Nt)

        Z_mod = Z_mod[I0:]
        T = T[I0:]
        Nt = len(T)

        var = np.std(Z_mod[:, JL])
        if var < 0.05:

            VR_mean = 0
            VL_mean = 0
            AR_mean = np.amax(Z_mod[:, JR])
            AL_mean = np.amax(Z_mod[:, JL])
            TR_mean = 0
            TL_mean = 0
        else:
            T_measure = 2
            VR_mean, AR_mean, TR_mean = mean_velocity(Z_mod, X, T, T_measure, JR)
            VL_mean, AL_mean, TL_mean = mean_velocity(Z_mod, X, T, T_measure, JL)
        velocities_R.append(VR_mean)
        velocities_L.append(VL_mean)
        amplitudes_R.append(AR_mean)
        amplitudes_L.append(AL_mean)
        periods_R.append(TR_mean)
        periods_L.append(TL_mean)
        gammas.append(gamma)
    plt.scatter(gammas, np.abs(velocities_R), c="b")
    plt.scatter(gammas, np.abs(velocities_L), c="r")
    plt.show()
    plt.close()

    plt.scatter(gammas, np.abs(amplitudes_R), c="b")
    plt.scatter(gammas, np.abs(amplitudes_L), c="r")
    plt.show()
    plt.close()

    plt.scatter(gammas, np.abs(periods_R), c="b")
    plt.scatter(gammas, np.abs(periods_R), c="r")
    plt.show()
    plt.close()
    np.savetxt(working_directory + '/parameter.txt', gammas, delimiter=',')
    np.savetxt(working_directory + '/v_R.txt', velocities_R, delimiter=',')
    np.savetxt(working_directory + '/v_L.txt', velocities_L, delimiter=',')
    np.savetxt(working_directory + '/A_R.txt', amplitudes_R, delimiter=',')
    np.savetxt(working_directory + '/A_L.txt', amplitudes_L, delimiter=',')
    np.savetxt(working_directory + '/T_R.txt', periods_R, delimiter=',')
    np.savetxt(working_directory + '/T_L.txt', periods_L, delimiter=',')