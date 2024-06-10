import matplotlib.pyplot as plt

from back_process import *

if __name__ == '__main__':
    directory = 'C:/mnustes_science/simulation_data/FD/PT_extended/variando_distancia_rabi/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.0180/sigma=6.000/gamma=0.1850/dist=24.0000'

    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    X = np.loadtxt(directory + '/X.txt', delimiter=',')
    T = np.loadtxt(directory + '/T.txt', delimiter=',')
    params = np.loadtxt(directory + '/parameters.txt', delimiter=',')
    eig_r = np.loadtxt(directory + '/eigenvalues_r.txt', delimiter=',')
    eig_i = np.loadtxt(directory + '/eigenvalues_i.txt', delimiter=',')
    # [alpha, beta, gamma_0, dist,  mu, nu, sigma]
    dist = params[3]

    Nt = len(T)
    each_T = 10
    ti = int(0.9 * Nt)

    Z_r = Z_r[ti::each_T, :]
    Z_i = Z_i[ti::each_T, :]
    T = T[ti::each_T]
    Nt = len(T)
    R_01 = []
    I_01 = []
    R_02 = []
    I_02 = []
    R_03 = []
    I_03 = []
    R_04 = []
    I_04 = []
    R_05 = []
    I_05 = []
    for i in range(Nt):
        x = np.argmax(eig_r[i, :])
        y = np.argsort(eig_r[i, :])
        r = np.take_along_axis(eig_r[i, :], y, axis=0)
        i = np.take_along_axis(eig_i[i, :], y, axis=0)
        R_01.append(r[-1])
        I_01.append(i[-1])
        R_02.append(r[-2])
        I_02.append(i[-2])
        R_03.append(r[-3])
        I_03.append(i[-3])
        R_04.append(r[-4])
        I_04.append(i[-4])
        R_05.append(r[-5])
        I_05.append(i[-5])

    def animate(i):
        line_1.set_offsets((R_01[i], I_01[i]))
        line_2.set_offsets((R_02[i], I_02[i]))
        line_3.set_offsets((R_03[i], I_03[i]))
        line_4.set_offsets((R_04[i], I_04[i]))
        line_5.set_offsets((R_05[i], I_05[i]))
        curve_1.set_data(X, Z_r[i, :])
        curve_2.set_data(X, Z_i[i, :])
        curve_mod.set_data(X, np.sqrt(Z_r[i, :] ** 2 + Z_i[i, :] ** 2))
        curve_1.set_color("b")
        curve_2.set_color('r')
        curve_mod.set_color("k")
        if i == 1:
            ax2.legend()
        return line_1, line_2, line_3, line_4, curve_1, curve_2, curve_mod

    fig, ((ax1, ax2)) = plt.subplots(nrows=2, ncols=1)

    line_1 = ax1.scatter(R_01[0], I_01[0], c="k")
    line_2 = ax1.scatter(R_02[0], I_02[0], c="k")
    line_3 = ax1.scatter(R_03[0], I_03[0], c="k")
    line_4 = ax1.scatter(R_04[0], I_04[0], c="k")
    line_5 = ax1.scatter(R_05[0], I_05[0], c="k")
    ax1.set_xlim(-0.1, 0.1)
    ax1.set_ylim(-0.2, 0.2)

    ax2.set_xlim(X[0], X[-1])
    ax2.set_ylim(-0.5, 0.5)

    curve_1, = ax2.plot(X, Z_r[0, :], lw=2, c="b", label="$A_R$")
    curve_2, = ax2.plot(X, Z_i[0, :], lw=2, c="r", label="$A_I$")
    curve_mod, = ax2.plot(X, np.sqrt(Z_i[0, :] ** 2 + Z_r[0, :] ** 2), lw=2, c="r", ls="dashed")
    ax1.hlines(0, -0.1, 0.1, colors="k", alpha=0.8)
    ax1.vlines(0, -0.2, 0.2, colors="k", alpha=0.8)

    ani = FuncAnimation(fig, animate,
                        frames=Nt,
                        interval=1,
                        blit=True)
    FFwriter = animation.FFMpegWriter()
    ani.save("test_02.gif", dpi=300, writer=PillowWriter(fps=120))
    plt.close()