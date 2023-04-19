from functions import *
from back_process import *
from time_integrators import *
from directories_lyap import *

if __name__ == '__main__':
    disco = 'E:/'
    initial_dir_data = str(disco) + "mnustes_science/simulation_data/FD/PDNLS_chaosa" + subdir_branching
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    H = np.loadtxt(directory + '/hamiltonian_t.txt', delimiter=',')
    T = np.loadtxt(directory + '/T.txt', delimiter=',')[0:-1]
    H_mean = filtro_array(300, H)#filtro_array(50, H)
    plt.plot(T, H)
    plt.plot(T, H_mean)
    plt.show()
    plt.close()
    print(len(T))
    print(len(H_mean))
    Nt = len(T)
    dt = T[1] - T[0]

    D = sparse_D_neumann(Nt, dt)
    DD = sparse_DD_neumann(Nt, dt)
    print(D.shape, H_mean.shape)
    D1_H = D.dot(np.transpose(H_mean))
    D2_H = DD.dot(np.transpose(H_mean))

    T_points_max = []
    H_points_max = []
    T_points_min = []
    H_points_min = []
    T_points = []
    H_points = []
    for i in range(Nt):
        if np.sign(D1_H[i]) != np.sign(D1_H[i - 1]) and D2_H[i] < 0 and i != 0:
            T_points_max.append(T[i])
            H_points_max.append(H_mean[i])
            T_points.append(T[i])
            H_points.append(H_mean[i])
        elif np.sign(D1_H[i]) != np.sign(D1_H[i - 1]) and D2_H[i] > 0 and i != 0:
            T_points_min.append(T[i])
            H_points_min.append(H_mean[i])
            T_points.append(T[i])
            H_points.append(H_mean[i])
    H_points_max = np.array(H_points_max)
    H_points_min = np.array(H_points_min)
    T_points_max = np.array(T_points_max)
    T_points_min = np.array(T_points_min)
    H_points = np.array(H_points)
    T_points = np.array(T_points)
    np.savetxt(directory + '/H_mean.txt', H_mean, delimiter=',')
    np.savetxt(directory + '/H_max.txt', H_points_max, delimiter=',')
    np.savetxt(directory + '/H_min.txt', H_points_min, delimiter=',')
    np.savetxt(directory + '/T_max.txt', T_points_max, delimiter=',')
    np.savetxt(directory + '/T_min.txt', T_points_min, delimiter=',')
    np.savetxt(directory + '/H_ext.txt', H_points, delimiter=',')
    np.savetxt(directory + '/T_ext.txt', T_points, delimiter=',')

    plt.scatter(T_points_max, H_points_max, color="r",zorder=1)
    plt.scatter(T_points_min, H_points_min, color="b", zorder=1)
    plt.plot(T, H_mean, color="k", zorder=0)
    #plt.plot(T, D1_H)
    #plt.plot(T, D2_H)
    plt.xlabel('$t$', size='20')
    plt.xticks(fontsize=15)
    plt.xlim([T[0], T[-1]])
    plt.ylabel('$H(t)$', size='20')
    plt.yticks(fontsize=15)
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.show()
    plt.close()