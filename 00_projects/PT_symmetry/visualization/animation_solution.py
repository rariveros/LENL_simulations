from back_process import *


if __name__ == '__main__':
    directory = 'C:/mnustes_science/simulation_data/FD/PT_extended/variando_distancia_rabi/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.0180/sigma=6.000/gamma=0.1850/dist=21.0000'

    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    X = np.loadtxt(directory + '/X.txt', delimiter=',')
    T = np.loadtxt(directory + '/T.txt', delimiter=',')
    params = np.loadtxt(directory + '/parameters.txt', delimiter=',')
    # [alpha, beta, gamma_0, dist,  mu, nu, sigma]
    dist = params[3]
    beta = 0.004811649356064012

    Z_r = Z_r / np.sqrt(beta)
    Z_i = Z_i / np.sqrt(beta)

    Nt = len(T)
    each_T = 40
    ti = int(0.3 * Nt)

    Z_r = Z_r[ti::each_T, :]
    Z_i = Z_i[ti::each_T, :]
    T = T[ti::each_T]

    Nt = len(T)

    def animate(i):
        line_1.set_data(X, Z_r[i, :])
        line_2.set_data(X, Z_i[i, :])
        line_3.set_data(X,  (Z_r[i, :] ** 2 + Z_i[i, :] ** 2) ** 0.5)
        return line_1, line_2, line_3

    def init():
        line_1.set_data([], [])
        line_2.set_data([], [])
        line_3.set_data([], [])
        line_1.set_color("b")
        line_2.set_color('r')
        line_3.set_color('k')
        return line_1, line_2, line_3

    fig = plt.figure()
    axis = plt.axes(xlim=(X[0], X[-1]),
                    ylim=(-7, 7))
    plt.xlabel('$x$', fontsize=20)
    plt.ylabel('$A_R, A_I$', fontsize=20)
    plt.grid(alpha=0.2, color="k")
    #plt.axis("off")
    axis.set_aspect(60 * np.sqrt(beta))

    line_1, = axis.plot([], [], lw=2)
    line_2, = axis.plot([], [], lw=2)
    line_3, = axis.plot([], [], lw=2, ls="--")

    ani = FuncAnimation(fig, animate,
                                   init_func=init,
                                   frames=Nt,
                                   interval=1,
                                   blit=True)
    FFwriter = animation.FFMpegWriter()
    ani.save("rabi_example.gif", dpi=300, writer=PillowWriter(fps=120))
    plt.close()

    #FFwriter = animation.FFMpegWriter()
    #anim.save('name_test.mp4', writer=FFwriter, dpi=300)
    #plt.close()