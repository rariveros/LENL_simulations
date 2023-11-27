from functions import *
from back_process import *
from time_integrators import *
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

if __name__ == '__main__':
    disco = 'E:/'
    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FD'
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    directories = [name for name in os.listdir(working_directory) if os.path.isdir(os.path.join(working_directory, name))]
    yf = []
    xf = []
    n = 0
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for i in range(len(directories)):
        directory = working_directory + "/" + directories[i]
        Z_r_points_np = np.loadtxt(directory + '/Z_r_points.txt', delimiter=',')
        x_points_np = np.loadtxt(directory + '/x_points.txt', delimiter=',')
        t_points_np = np.loadtxt(directory + '/t_points.txt', delimiter=',')
        n = n + 1
        Nt = len(t_points_np)
        dt = np.abs(t_points_np[1] - t_points_np[0])
        yf_i = fft(Z_r_points_np)
        yf_i = 2.0 / Nt * np.abs(yf_i[0:Nt // 2])
        xf_i = fftfreq(Nt, dt)[:Nt // 2]
        yf_i = yf_i[1:]
        xf_i = xf_i[1:]
        yf.append(yf_i)
        xf.append(xf_i)
        #plt.plot(xf_i, np.log(yf_i), color="k", linewidth=0.5)
        #plt.xlim(0, 2)
        #plt.xlabel("$\\textrm{Frequency }(\Omega)$", fontsize=28)
        #plt.ylabel("$\\textrm{Power (dB)}$", fontsize=28)
        #plt.xticks(fontsize=15)
        #plt.yticks(fontsize=15)
        #plt.grid(alpha=0.5)
        #plt.tight_layout()
        #plt.savefig(subdirectory + '/fft_amplitude.png', dpi=300)
        #plt.close()
        gamma_i = float(directories[i].split("=")[-1])
        color = [(gamma_i - 0.15) / (0.2 - 0.15), 0, 0]
        alpha = (2 - (gamma_i - 0.15) / (0.2 - 0.15)) / 2
        ax.plot(gamma_i * np.ones(len(xf_i[0:int(0.2*len(xf_i))])), xf_i[0:int(0.2*len(xf_i))], np.log(yf_i)[0:int(0.2*len(xf_i))], color=color, alpha=alpha)
    ax.set_xlabel('$\sigma_i$', fontsize=22)
    ax.set_ylabel('$\gamma_0$', fontsize=22)
    ax.zaxis.set_rotate_label(False)
    ax.set_zlabel('$P^{\\frac{1}{2}}$', fontsize=22, rotation=0)
    ax.axes.set_xlim3d(left=0.15, right=0.2)
    ax.axes.set_ylim3d(bottom=0, top=1)
    ax.axes.set_zlim3d(bottom=-13)
    ax.axes.tick_params(axis="both", labelsize=12)
    ax.view_init(elev=21, azim=-145, roll=0)
    plt.show()
