from back_process import *

if __name__ == '__main__':

    #fig, ax = plt.subplots()
    #ax.set_xlim([-5, 5])
    #ax.set_ylim([-5, 0])
    t_grid = np.arange(0, 1, 0.001)
    freq1 = 20
    freq2 = 1
    x0 = 0.5
    x10 = -1
    x20 = 1
    R = 2
    #scat_1 = ax.scatter(1, 0, s=300, c=np.array([68, 50, 230, 255]) / 255, zorder=10)
    #scat_1 = ax.scatter(1, 0, s=300, c=np.array([68, 50, 230, 255]) / 255, zorder=10)
    #scat_2 = ax.scatter(1, 0, s=300, c=np.array([68, 50, 230, 255]) / 255, zorder=10)
    #line_1, = ax.plot([], [], lw=3, color="k", zorder=2)
    #line_2, = ax.plot([], [], lw=3, color="k", zorder=2)
    x1 = x0 * np.cos(2 * np.pi * freq1 * t_grid) * np.cos(2 * np.pi * freq2 * t_grid + np.pi / 2)
    y1 = -np.sqrt(R ** 2 - x1 ** 2)
    x2 = x0 * np.cos(2 * np.pi * freq1 * t_grid) * np.cos(2 * np.pi * freq2 * t_grid)
    y2 = -np.sqrt(R ** 2 - x2 ** 2)


    #def animate(i):
    #    scat_1.set_offsets((x1[i] + x10, y1[i]))
    #    scat_2.set_offsets((x2[i] + x20, y2[i]))
    #    line_1.set_data(np.array([x10, x1[i] + x10]), np.array([0, y1[i]]))
    #    line_2.set_data(np.array([x20, x2[i] + x20]), np.array([0, y2[i]]))
    #    return scat_1, scat_2, line_1, line_2,


    #ani = animation.FuncAnimation(fig, animate, repeat=True,
    #                              frames=len(x1) - 1, interval=1)
    #ani.save("unbroken_PT.gif", dpi=300, writer=PillowWriter(fps=60))
    #plt.close()

    # plot 1:
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
    ax1.plot(t_grid, x1, color=np.array([230, 50, 50, 255]) / 255)
    ax2.plot(t_grid, x2, color=np.array([68, 50, 230, 255]) / 255)
    ax1.set_xlim([t_grid[0], t_grid[-1]])
    ax2.set_xlim([t_grid[0], t_grid[-1]])

    plt.savefig("test_01.png", dpi=300)


    #FFwriter = animation.FFMpegWriter()
    #anim.save('name_test.mp4', writer=FFwriter, dpi=300)
    #plt.close()