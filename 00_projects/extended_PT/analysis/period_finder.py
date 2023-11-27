import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    disco = 'E:/'

    initial_dir_data = str(disco) + "mnustes_science/simulation_data/FD/PDNLS_chaosa"
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    analysis_directory = directory + "//analysis"
    S_plus = np.loadtxt(analysis_directory + '/S+.txt', delimiter=',')
    S_minus = np.loadtxt(analysis_directory + '/S-.txt', delimiter=',')
    T = np.loadtxt(analysis_directory + '/tau.txt', delimiter=',')
    X = np.loadtxt(directory + '/X.txt', delimiter=',')
    #params = np.loadtxt(directory + '/params.txt', delimiter=',')
    #gaussian_L = params[1]
    #gaussian_R = params[2]

    gaussian_L_j = np.where(X == -37)[0][0]
    gaussian_R_j = np.where(X == 37)[0][0]
    S_L = filtro_array(10, S_plus[:, gaussian_L_j])
    S_R = filtro_array(10, S_plus[:, gaussian_R_j])

    #S_L = S_plus[:, gaussian_L_j]
    #S_R = S_plus[:, gaussian_R_j]

    Nt = len(T)
    dt = T[1] - T[0]

    D = sparse_D_neumann(Nt, dt)
    DD = sparse_DD_neumann(Nt, dt)
    D1_SL = D.dot(np.transpose(S_L))
    D2_SL = DD.dot(np.transpose(S_L))
    D1_SR = D.dot(np.transpose(S_R))
    D2_SR = DD.dot(np.transpose(S_R))

    T_L_points_max = []
    S_L_points_max = []
    T_R_points_max = []
    S_R_points_max = []
    SL_max = np.amax(S_L)
    SR_max = np.amax(S_R)
    for i in range(Nt):
        if np.sign(D1_SL[i]) != np.sign(D1_SL[i - 1]) and D2_SL[i] < 0 and i != 0:
            T_L_points_max.append(T[i])
            S_L_points_max.append(S_L[i])
        elif np.sign(D1_SR[i]) != np.sign(D1_SR[i - 1]) and D2_SR[i] < 0 and i != 0:
            T_R_points_max.append(T[i])
            S_R_points_max.append(S_R[i])
    S_L_points_max = np.array(S_L_points_max)
    T_L_points_max = np.array(T_L_points_max)
    S_R_points_max = np.array(S_R_points_max)
    T_R_points_max = np.array(T_R_points_max)
    dT_L = np.diff(T_L_points_max[1:-1])
    dT_R = np.diff(T_R_points_max[1:-1])

    period_L = np.mean(dT_L)
    period_L_std = np.std(dT_L)
    period_R = np.mean(dT_R)
    period_R_std = np.std(dT_R)
    print(1/period_L)
    print(period_L_std)
    print(1/period_R)
    print(period_R_std)

    np.savetxt(analysis_directory + '/S_L.txt', S_L, delimiter=',')
    np.savetxt(analysis_directory + '/S_L_points_max.txt', S_L_points_max, delimiter=',')
    np.savetxt(analysis_directory + '/T_L_max.txt', T_L_points_max, delimiter=',')
    np.savetxt(analysis_directory + '/S_R.txt', S_R, delimiter=',')
    np.savetxt(analysis_directory + '/S_R_points_max.txt', S_R_points_max, delimiter=',')
    np.savetxt(analysis_directory + '/T_R_max.txt', T_R_points_max, delimiter=',')

    plt.vlines(T_L_points_max, 0, 1.1 * np.amax(S_L), color="r", zorder=1)
    plt.plot(T, S_plus[:, gaussian_L_j], color="k", zorder=0)
    plt.xlabel('$t$', size='20')
    plt.xticks(fontsize=15)
    plt.xlim([T[0], T[-1]])
    plt.ylabel('$S_{+}(x_L, t)$', size='20')
    plt.yticks(fontsize=15)
    plt.ylim([0, 1.1 * np.amax(S_L)])
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.show()
    #plt.savefig(analysis_directory + '/S+_L_periods.png', dpi=300)
    plt.close()

    plt.vlines(T_R_points_max, 0, 1.1 * np.amax(S_R), color="r", zorder=1)
    plt.plot(T, S_plus[:, gaussian_R_j], color="k", zorder=0)
    plt.xlabel('$t$', size='20')
    plt.xticks(fontsize=15)
    plt.xlim([T[0], T[-1]])
    plt.ylabel('$S_{+}(x_R, t)$', size='20')
    plt.yticks(fontsize=15)
    plt.ylim([0, 1.1 * np.amax(S_R)])
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.show()
    #plt.savefig(analysis_directory + '/S+_R_periods.png', dpi=300)
    plt.close()