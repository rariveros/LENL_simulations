from back_process import *

if __name__ == '__main__':
    Lx = 120
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

    X = np.arange(-Lx / 2, Lx / 2, 0.05)


    Z_1 = a * np.cos(k * X)
    Z_2 = a * np.ones(len(X))

    fig = plt.figure()
    axis = plt.axes(xlim=(X[0], X[-1]),
                    ylim=(-1.5, 1.5))
    #plt.xlabel('$x$', fontsize=20)
    #plt.ylabel('$A_R(x)$', fontsize=20)
    plt.grid(alpha=0.4)
    plt.plot(X, Z_1, c=(215 / 255, 2 / 255, 28 / 255), lw=2)
    plt.plot(X, Z_2, c='k', lw=2, ls="dashed")
    #plt.axis("off")
    axis.set_aspect(2)

    plt.savefig('C_example_homo.png', dpi=300)
    plt.close()