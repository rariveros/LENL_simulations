from back_process import *

if __name__ == '__main__':
    disco = 'D:/'
    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FD/PDNLS_oscillatory/alpha=5.721/beta=1.000/mu=0.100/nu=0.056/sigma=60.000'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    X = np.loadtxt(directory + '/X.txt', delimiter=',')
    T = np.loadtxt(directory + '/T.txt', delimiter=',')
    Nt = len(T)

    t_i = int(0.9 * Nt)
    t_f = -1

    Z_r = Z_r[t_i:t_f, :] / np.sqrt(0.005)
    Z_i = Z_i[t_i:t_f, :] / np.sqrt(0.005)
    T = T[t_i: t_f]
    Nt = len(T)

    def animate(i):
        line_1.set_data(X, Z_r[i, :])
        #line_2.set_data(X, Z_i[i, :])
        return line_1,#, line_2

    def init():
        line_1.set_data([], [])
        #line_2.set_data([], [])
        line_1.set_color('k')
        #line_2.set_color((215 / 255, 2 / 255, 28 / 255))
        return line_1,#, line_2

    fig = plt.figure()
    axis = plt.axes(xlim=(X[0], X[-1]),
                    ylim=(1.1 * np.amin(Z_r), 1.1 * np.amax(Z_r)))
    plt.xlabel('$x\ (\\textrm{mm})$', fontsize=15)
    plt.ylabel('$A_R(x)$', fontsize=15)
    plt.grid(alpha=0.4)
    #plt.axis("off")
    axis.set_aspect(1)

    line_1, = axis.plot([], [], lw=1)
    #line_2, = axis.plot([], [], lw=1)

    ani = FuncAnimation(fig, animate,
                                   init_func=init,
                                   frames=Nt,
                                   interval=1,
                                   blit=True)
    FFwriter = animation.FFMpegWriter(fps=60)
    ani.save(directory + '/real_field_animation.mp4', writer=FFwriter, dpi=300)
    plt.close()

    #FFwriter = animation.FFMpegWriter()
    #anim.save('name_test.mp4', writer=FFwriter, dpi=300)
    #plt.close()