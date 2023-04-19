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
    sigma = 10
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

    X_short = X[sigma_window_l:sigma_window_r]
    Z_modulo_short = Z_modulo[:, sigma_window_l:sigma_window_r]
    Z_r_short = Z_r[:, sigma_window_l:sigma_window_r]
    Z_i_short = Z_i[:, sigma_window_l:sigma_window_r]
    D1_Z_short = D1_Z[:, sigma_window_l: sigma_window_r]
    D2_Z_short = D2_Z[:, sigma_window_l: sigma_window_r]
    Nx = len(Z_modulo[0, :])

    Z_r_points = []
    Z_i_points = []
    amplitude_points = []
    x_points = []
    t_points = []
    for i in range(Nt):
        x_points_j = []
        amplitude_points_j = []
        t_points_i = []
        Z_r_points_j = []
        Z_i_points_j = []
        for j in range(sigma_window_l, sigma_window_r):
            if np.sign(D1_Z[i, j]) != np.sign(D1_Z[i, j - 1]) and D2_Z[i, j] < 0 and j != 0:
                Z_r_points_j.append(Z_r[i, j])
                Z_i_points_j.append(Z_i[i, j])
                amplitude_points_j.append(Z_modulo[i, j])
                x_points_j.append(X[j])
                t_points_i.append(T[i])
        Z_r_points.append(Z_r_points_j)
        Z_i_points.append(Z_i_points_j)
        amplitude_points.append(amplitude_points_j)
        t_points.append(t_points_i)
        x_points.append(x_points_j)


    x_points = np.array(x_points)
    t_points = np.array(t_points)
    amplitude_points = np.array(amplitude_points)
    amplitude_points = amplitude_points - np.mean(amplitude_points)
    Z_r_points = np.array(Z_r_points)
    Z_r_points = Z_r_points - np.mean(Z_r_points)
    Z_i_points = np.array(Z_i_points)
    Z_i_points = Z_i_points - np.mean(Z_i_points)

    Z_modules_FFT = []
    for j in range(len(amplitude_points[0, :])):
        Z_module_FFT = fft(amplitude_points[:, j])
        Z_module_FFT = 2.0 / Nt * np.abs(Z_module_FFT[0:Nt // 2])
        xf_i = fftfreq(Nt, dt)[:Nt // 2]
        Z_module_FFT = np.array(Z_module_FFT)
        Z_modules_FFT.append(Z_module_FFT)
    Z_modules_FFT = np.transpose(np.array(Z_modules_FFT))

    fig, ax = plt.subplots()

    # legend
    pcm = plt.pcolormesh(X, T, Z_modulo, cmap=parula_map, vmin=np.amin(Z_modulo), vmax=np.amax(Z_modulo), shading='auto')
    cbar = plt.colorbar(pcm)
    cbar.set_label('$|A|$', rotation=0, size=25, labelpad=-27, y=1.11)
    print(x_points[:, 0])
    print(t_points[:, 0])
    ax.plot(x_points, t_points, zorder=1, linewidth=1, color="r")
    # put the major ticks at the middle of each cell
    plt.xlim(-10, 10)
    plt.xlabel('$x$', size='25')
    plt.xticks(fontsize=15)

    plt.ylim([1800, T[-1]])
    plt.ylabel('$t$', size='25')
    plt.yticks(fontsize=15)

    plt.grid(linestyle='--', alpha=0.2, color='k')
    cbar.ax.tick_params(labelsize=15)

    # labels
    plt.tight_layout()
    plt.savefig(directory + '/field_visual_zoom.png', dpi=300)
    plt.close()

    fig, (ax1, ax2, ax3) = plt.subplots(3)

    ax1.plot(T, amplitude_points[:, 0] - np.mean(amplitude_points[:, 0]))
    ax1.set_ylabel("$\delta |A(x_1, t)|$", fontsize="12")
    ax1.set_xlim(T[-400], T[-1])
    ax2.plot(T, amplitude_points[:, 1] - np.mean(amplitude_points[:, 1]), color="green")
    ax2.set_ylabel("$\delta |A(x_2, t)|$", fontsize="12")
    ax2.set_xlim(T[-400], T[-1])
    ax3.plot(T, amplitude_points[:, 2] - np.mean(amplitude_points[:, 2]), color="red")
    ax3.set_ylabel("$\delta |A(x_3, t)|$", fontsize="12")
    ax3.set_xlim(T[-400], T[-1])
    #ax4.plot(T, amplitude_points[:, 3] - np.mean(amplitude_points[:, 2]), color="c")
    #ax4.set_ylabel("$\delta |A(x_4, t)|$", fontsize="12")
    #ax4.set_xlim(T[-400], T[-1])
    #ax5.plot(T, amplitude_points[:, 4] - np.mean(amplitude_points[:, 2]), color="m")
    #ax5.set_ylabel("$\delta |A(x_5, t)|$", fontsize="12")
    #ax5.set_xlim(T[-400], T[-1])
    #ax5.set_xlabel("$t$", fontsize="15")
    ax1.grid(alpha=0.2)
    ax2.grid(alpha=0.2)
    ax3.grid(alpha=0.2)
    #ax4.grid(alpha=0.2)
    #ax5.grid(alpha=0.2)
    plt.sca(ax1)
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.sca(ax2)
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    #plt.sca(ax3)
    #plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    #plt.sca(ax4)
    #plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.savefig(directory + '/oscillations.png', dpi=300)
    plt.close()

    fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5)
    maximums = int(len(amplitude_points[0, :]) / 2)
    ax1.plot(xf_i, Z_modules_FFT[:, 0])
    ax1.set_ylabel("$\delta |A(x_1, \omega)|$", fontsize="12")
    ax1.set_xlim(0, 1)
    ax2.plot(xf_i, Z_modules_FFT[:, 1], color="green")
    ax2.set_ylabel("$\delta |A(x_2, \omega)|$", fontsize="12")
    ax2.set_xlim(0, 1)
    ax3.plot(xf_i, Z_modules_FFT[:, 2], color="red")
    ax3.set_ylabel("$\delta |A(x_3, \omega)|$", fontsize="12")
    ax3.set_xlim(0, 1)
    #ax4.plot(xf_i, Z_modules_FFT[:, 3], color="c")
    #ax4.set_ylabel("$\delta |A(x_4, \omega)|$", fontsize="12")
    #ax4.set_xlim(0, 1)
    #ax5.plot(xf_i, Z_modules_FFT[:, 4], color="m")
    #ax5.set_ylabel("$\delta |A(x_5, \omega)|$", fontsize="12")
    #ax5.set_xlim(0, 1)
    #ax5.set_xlabel("$\omega$", fontsize="15")
    ax1.grid(alpha=0.2)
    ax2.grid(alpha=0.2)
    ax3.grid(alpha=0.2)
    #ax4.grid(alpha=0.2)
    #ax5.grid(alpha=0.2)
    plt.sca(ax1)
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.sca(ax2)
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    #plt.sca(ax3)
    #plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    #plt.sca(ax4)
    #plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.savefig(directory + '/FFTs.png', dpi=300)
    plt.close()

    plt.scatter(Z_i_points[:, 0], Z_r_points[:, 0], s=1.5)
    plt.xlabel("$\delta A_R(x_1, t)$", fontsize=20)
    plt.ylabel("$\delta A_I(x_1, t)$", fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig(directory + '/phase_01.png', dpi=300)
    plt.close()

    plt.scatter(Z_i_points[:, 1], Z_r_points[:, 1], s=1.5, color="green")
    plt.xlabel("$\delta A_R(x_2, t)$", fontsize=20)
    plt.ylabel("$\delta A_I(x_2, t)$", fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig(directory + '/phase_02.png', dpi=300)
    plt.close()

    plt.scatter(Z_i_points[:, 2], Z_r_points[:, 2], s=1.5, color="red")
    plt.xlabel("$\delta A_R(x_3, t)$", fontsize=20)
    plt.ylabel("$\delta A_I(x_3, t)$", fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig(directory + '/phase_03.png', dpi=300)
    plt.close()