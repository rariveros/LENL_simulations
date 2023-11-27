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

    T_points_max = []
    H_points_max = []
    T_points_min = []
    H_points_min = []
    T_points = []
    H_points = []
    n = 0
    files = [name for name in os.listdir(directory) if os.path.isfile(os.path.join(directory, name))]
    for file in files:
        if file[:3] == "ham":
            H = np.loadtxt(directory + '/' + file, delimiter=',')
            T = np.loadtxt(directory + '/T.txt', delimiter=',')
            H_mean = filtro_array(100, H)#filtro_array(50, H)
            Nt = len(T)
            dt = T[1] - T[0]

            D = sparse_D_neumann(Nt, dt)
            DD = sparse_DD_neumann(Nt, dt)
            print(D.shape, H_mean.shape)
            D1_H = D.dot(np.transpose(H_mean))
            D2_H = DD.dot(np.transpose(H_mean))

            for i in range(Nt):
                n = n + 1
                if np.sign(D1_H[i]) != np.sign(D1_H[i - 1]) and D2_H[i] < 0 and i != 0:
                    T_points_max.append(n * dt)
                    H_points_max.append(H_mean[i])
                    T_points.append(n * dt)
                    H_points.append(H_mean[i])
                elif np.sign(D1_H[i]) != np.sign(D1_H[i - 1]) and D2_H[i] > 0 and i != 0:
                    T_points_min.append(n * dt)
                    H_points_min.append(H_mean[i])
                    T_points.append(n * dt)
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