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
    dists = []
    lamb_L = []
    lamb_R = []
    lamb_L_std = []
    lamb_R_std = []
    for directory in directories:
        Z_r = np.loadtxt(working_directory + "/" + directory + '/field_real.txt', delimiter=',')
        Z_i = np.loadtxt(working_directory + "/" + directory + '/field_img.txt', delimiter=',')
        X = np.loadtxt(working_directory + "/" + directory + '/X.txt', delimiter=',')
        T = np.loadtxt(working_directory + "/" + directory + '/T.txt', delimiter=',')
        params = np.loadtxt(working_directory + "/" + directory + '/parameters.txt', delimiter=',')
        print("######### " + directory + " #########")
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
        Z_mod = np.abs(Z)
        lamb_L_i = []
        lamb_R_i = []
        lamb_L_std_i = []
        lamb_R_std_i = []
        for j in range(Nt):
            Z_mod_L = Z_mod[j, :int(Nx / 2)]
            Z_mod_R = Z_mod[j, int(Nx / 2):]
            X_L = X[:int(Nx / 2)]
            X_R = X[int(Nx / 2):]
            arg_max_L = argrelextrema(Z_mod_L, np.greater)
            arg_max_R = argrelextrema(Z_mod_R, np.greater)
            lambdas_L_j = 2 * np.diff(X_L[arg_max_L])
            lambdas_R_j = 2 * np.diff(X_R[arg_max_R])
            lamb_L_i.append(np.mean(lambdas_L_j))
            lamb_R_i.append(np.mean(lambdas_R_j))
        lamb_L.append(np.mean(lamb_L_i))
        lamb_R.append(np.mean(lamb_R_i))
        lamb_L_std.append(np.std(lamb_L_i))
        lamb_R_std.append(np.std(lamb_R_i))
        dists.append(params[3])
    np.savetxt('lamb_L.txt', lamb_L, delimiter=',')
    np.savetxt('lamb_R.txt', lamb_R, delimiter=',')
    np.savetxt('lamb_L_std.txt', lamb_L_std, delimiter=',')
    np.savetxt('lamb_R_std.txt', lamb_R_std, delimiter=',')
    np.savetxt('dist.txt', dists, delimiter=',')
    plt.errorbar(dists, lamb_L, lamb_L_std, marker='o', ls='', ecolor="k", mec='black', color="r", label="left")
    plt.errorbar(dists, lamb_R, lamb_R_std, marker='o', ls='', ecolor="k", mec='black', color="b", label="right")
    plt.legend(fontsize=15)
    plt.xlabel("$d$", fontsize=15)
    plt.ylabel("$\lambda$", fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.show()
    plt.close()