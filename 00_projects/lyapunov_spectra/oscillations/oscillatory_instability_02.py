from functions import *
from back_process import *
from time_integrators import *
from directories_lyap import *

if __name__ == '__main__':
    disco = 'E:/'
    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FD' + main_directory + directoy_sim + subdir_branching
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    X = np.loadtxt(directory + '/X.txt', delimiter=',')
    T = np.loadtxt(directory + '/T.txt', delimiter=',')

    Nx = len(X)
    Nt = len(T)
    dt = T[1] - T[0]
    dx = X[1] - X[0]
    L = X[-1] - X[0]
    sigma_i = 60
    mu = 0.1
    nu = 0.056
    alpha = 5.721
    sigma = sigma_i #2.3548 * sigma_i#np.sqrt(np.sqrt(nu * alpha) * (sigma_i / mu))
    print(sigma)
    sigma_window_l = int((Nx - (sigma / L) * Nx) / 2)
    sigma_window_r = int(sigma_window_l + (sigma / L) * Nx)

    # DEFINE REGION OF INTEREST +- SIGMA/2

    ti = int(0.5 * Nt)
    tf = -1
    T = T[ti:tf]
    Nt = len(T)
    Z_r = Z_r[ti:tf, :]
    Z_i = Z_i[ti:tf, :]

    Z_complex = Z_r + 1j * Z_i

    # Definiendo variables finales
    Z_modulo= np.absolute(Z_complex)
    Z_arg = np.angle(Z_complex)
    Z_arg = (2 * np.pi + Z_arg) * (Z_arg < 0) + Z_arg * (Z_arg > 0)

    D = sparse_D_neumann(Nx, dx)
    DD = sparse_DD_neumann(Nx, dx)

    D1_Z = Der(D, np.transpose(Z_modulo))
    D2_Z = Der(DD, np.transpose(Z_modulo))
    D1_Z = np.transpose(D1_Z)
    D2_Z = np.transpose(D2_Z)

    J_points_init = []
    K = 0
    for j in range(sigma_window_l, sigma_window_r):
        if np.sign(D1_Z[0, j]) != np.sign(D1_Z[0, j - 1]) and D2_Z[0, j] < 0 and j != 0:
            K = K + 1
            J_points_init.append(j)
    print(K)

    Z_r_points = []
    Z_i_points = []
    amplitude_points = []
    x_points = []
    t_points = []
    window_size = 4
    for k in range(K):
        window_left = J_points_init[k] - int(window_size / 2)
        window_right = J_points_init[k] + int(window_size / 2)
        x_points_k = []
        amplitude_points_k = []
        t_points_k = []
        Z_r_points_k = []
        Z_i_points_k = []
        for i in range(Nt):
            for j in range(window_size):
                if np.sign(D1_Z[i, window_left + j]) != np.sign(D1_Z[i, window_left + j - 1]) and D2_Z[i, window_left + j] < 0 and j != 0:
                    Z_r_points_k.append(Z_r[i, window_left + j])
                    Z_i_points_k.append(Z_i[i, window_left + j])
                    amplitude_points_k.append(Z_modulo[i, window_left + j])
                    x_points_k.append(X[window_left + j])
                    t_points_k.append(T[i])
                    window_left = window_left + j - int(window_size / 2)
                    window_right = window_left + j + int(window_size / 2)
        Z_r_points.append(Z_r_points_k)
        Z_i_points.append(Z_i_points_k)
        amplitude_points.append(amplitude_points_k)
        t_points.append(t_points_k)
        x_points.append(x_points_k)



    Z_r_points_np = np.array(Z_r_points[2][:])
    Z_i_points_np = np.array(Z_i_points[2][:])
    x_points_np = np.array(x_points[2][:])
    t_points_np = np.array(t_points[2][:])
    np.savetxt(directory + '/Z_r_points.txt', Z_r_points_np, delimiter=',')
    np.savetxt(directory + '/Z_i_points.txt', Z_i_points_np, delimiter=',')
    np.savetxt(directory + '/x_points.txt', x_points_np, delimiter=',')
    np.savetxt(directory + '/t_points.txt', t_points_np, delimiter=',')


    fig, ax = plt.subplots()

    # legend
    pcm = plt.pcolormesh(X, T, Z_r, cmap=parula_map, vmin=-np.amax(Z_r), vmax=np.amax(Z_r), shading='auto')
    cbar = plt.colorbar(pcm)
    cbar.set_label('$A_R(x,t)$', rotation=0, size=25, labelpad=-27, y=1.11)
    for k in range(K):
        ax.plot(x_points[k], t_points[k], zorder=1, linewidth=1, color="r")
    # put the major ticks at the middle of each cell
    plt.xlim(X[0], X[-1])
    plt.xlabel('$x$', size='25')
    plt.xticks(fontsize=15)

    plt.ylim([T[0], T[-1]])
    plt.ylabel('$t$', size='25')
    plt.yticks(fontsize=15)

    plt.grid(linestyle='--', alpha=0.2, color='k')
    cbar.ax.tick_params(labelsize=15)

    # labels
    plt.tight_layout()
    plt.savefig(directory + '/test_01.png', dpi=300)
    plt.close()