from back_process import *


if __name__ == '__main__':
    directory = 'C:/mnustes_science/simulation_data/FD/PT_extended/variando_distancia/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.0180/sigma=6.000/gamma=0.1850/dist=30.2930'

    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    X = np.loadtxt(directory + '/X.txt', delimiter=',')
    T = np.loadtxt(directory + '/T.txt', delimiter=',')
    params = np.loadtxt(directory + '/parameters.txt', delimiter=',')
    Z_mod_rabi = np.loadtxt("mod_strobo_rabi.txt", delimiter=',')
    # [alpha, beta, gamma_0, dist,  mu, nu, sigma]
    dist = params[3]
    beta = 0.004811649356064012

    Z_r = Z_r / np.sqrt(beta)
    Z_i = Z_i / np.sqrt(beta)

    Z = Z_r + 1j * Z_i
    Z_conj = Z_r - 1j * Z_i
    Z_modulo = np.absolute(Z)
    fig = plt.figure()
    axis = plt.axes(xlim=(X[0], X[-1]),
                    ylim=(-6, 6))
    plt.xlabel('$x$', fontsize=20)
    plt.ylabel('$A_R, A_I$', fontsize=20)
    plt.grid(alpha=0.2, color="k")
    #plt.axis("off")
    axis.set_aspect(90 * np.sqrt(beta))

    line_1, = axis.plot([], [], lw=2)
    line_2, = axis.plot([], [], lw=2)
    line_3, = axis.plot([], [], lw=2, ls="--")
    line_1.set_data(X, 0 * Z_r[-1, :])
    line_2.set_data(X, 0 * Z_i[-1, :])
    line_3.set_data(X, 0 * (Z_r[-1, :] ** 2 + Z_i[-1, :] ** 2) ** 0.5)
    line_1.set_color("b")
    line_2.set_color('r')
    line_3.set_color('k')

    plt.savefig("quiscent.png", dpi=300)
    plt.close()