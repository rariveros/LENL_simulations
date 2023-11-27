from back_process import *
from matplotlib.animation import PillowWriter

if __name__ == '__main__':
    gamma_r = 1
    gamma_l = - gamma_r
    phi = 0 * np.pi
    dist = 1 * 33
    sigma = 7.9
    w = 2

    x_grid = np.arange(-50, 50, 0.1)
    t_grid = np.arange(0, 2 * np.pi / w, 0.1)
    Nx = len(x_grid)
    Nt = len(t_grid)
    X, T = np.meshgrid(x_grid, t_grid)

    y_grid = np.cos(w * T) * (gamma_l * np.exp(- (X - dist / 2) ** 2 / (2 * sigma ** 2)) + gamma_r * np.cos(phi) * np.exp(
        - (X + dist / 2) ** 2 / (2 * sigma ** 2)))

    def animate(i):
        line_1.set_data(x_grid, y_grid[i, :])
        line_2.set_data(x_grid, y_grid[0, :])
        return line_1, line_2,

    def init():
        line_1.set_data(x_grid, y_grid[0, :])
        line_1.set_color("k")
        line_2.set_color("r")
        return line_1, line_2,

    fig = plt.figure()
    axis = plt.axes(xlim=(x_grid[0], x_grid[-1]),
                    ylim=(-1.5, 1.5))
    #plt.xlabel('$x$', fontsize=20)
    #plt.ylabel('$u(x)$', fontsize=20)
    #plt.grid(alpha=0.4)
    #plt.axis("off")
    plt.tick_params(left=False, right=False, labelleft=False, labelbottom=False, bottom=False)
    plt.grid(alpha=0.3)
    axis.set_aspect(10)

    line_1, = axis.plot([], [], lw=3)
    line_2, = axis.plot([], [], lw=3)

    ani = FuncAnimation(fig, animate,
                                   init_func=init,
                                   frames=Nt,
                                   interval=1,
                                   blit=True)
    #FFwriter = animation.FFMpegWriter()
    #ani.save('name_test.gif', writer=FFwriter, dpi=300)
    writer = PillowWriter(fps=60)
    ani.save("demo2.gif", writer=writer)
    plt.close()

    #FFwriter = animation.FFMpegWriter()
    #anim.save('name_test.mp4', writer=FFwriter, dpi=300)
    #plt.close()