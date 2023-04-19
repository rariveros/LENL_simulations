from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    disco = 'E:/'
    initial_dir_data = str(disco) + "mnustes_science/simulation_data/FD/PDNLS_chaos/branching_dynamics"
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    Z_r = np.loadtxt(directory + '/field_real_FFT.txt', delimiter=',')
    K = np.loadtxt(directory + '/K.txt', delimiter=',')
    T = np.loadtxt(directory + '/T.txt', delimiter=',')
    Nx = len(K)
    Nt = len(T)
    ti = -100
    tf = -1
    Z_r = Z_r[ti:tf, :]
    T = T[ti:tf]
    Nt = len(T)

    def animate(i):
        line_1.set_data(K, Z_r[i, :])
        return line_1,

    def init():
        line_1.set_data([], [])
        line_1.set_color('k')
        return line_1,

    fig = plt.figure()
    axis = plt.axes(xlim=(0, 0.5),
                    ylim=(0, 0.5))
    plt.xlabel('$x$', fontsize=20)
    plt.ylabel('$\hat{A}(k)$', fontsize=20)
    plt.grid(alpha=0.4)

    line_1, = axis.plot([], [], lw=3)
    line_2, = axis.plot([], [], lw=3, ls="dashed")

    ani = FuncAnimation(fig, animate,
                                   init_func=init,
                                   frames=Nt,
                                   interval=0.01,
                                   blit=True)
    FFwriter = animation.FFMpegWriter()
    ani.save(directory + '/real_FFT.mp4', writer=FFwriter, dpi=300)
    plt.close()