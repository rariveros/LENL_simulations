from back_process import *


if __name__ == '__main__':
    Lx = 20
    x = np.arange(-Lx / 2, Lx / 2, 0.05)
    w = 4 / 3
    t = np.arange(0, 2 * np.pi / w, 0.1)
    Nx = len(x)
    Nt = len(t)

    a = 1
    b = 0.2
    n = 20
    k = 4
    sigma = 2

    X, T = np.meshgrid(x, t)
    Z_1 = a * np.cos(w * T + k * X) * np.exp(- X ** 2 / (2 * sigma) ** 2)
    Z_2 = a * np.exp(- X ** 2 / (2 * sigma) ** 2)

    def animate(i):
        line_1.set_data(x, Z_1[i, :])
        line_2.set_data(x, Z_2[i, :])
        return line_1, line_2

    def init():
        line_1.set_data([], [])
        line_2.set_data([], [])
        line_1.set_color((215 / 255, 2 / 255, 28 / 255))
        line_2.set_color('k')
        return line_1, line_2

    fig = plt.figure()
    axis = plt.axes(xlim=(x[0], x[-1]),
                    ylim=(-1.5, 1.5))
    #plt.xlabel('$x$', fontsize=20)
    #plt.ylabel('$u(x)$', fontsize=20)
    #plt.grid(alpha=0.4)
    plt.axis("off")
    axis.set_aspect(2)

    line_1, = axis.plot([], [], lw=3)
    line_2, = axis.plot([], [], lw=3, ls="dashed")

    ani = FuncAnimation(fig, animate,
                                   init_func=init,
                                   frames=Nt,
                                   interval=30,
                                   blit=True)
    FFwriter = animation.FFMpegWriter()
    ani.save("test_02.gif", dpi=300, writer=PillowWriter(fps=120))
    plt.close()

    #FFwriter = animation.FFMpegWriter()
    #anim.save('name_test.mp4', writer=FFwriter, dpi=300)
    #plt.close()