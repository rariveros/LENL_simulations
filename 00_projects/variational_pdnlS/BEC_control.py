from functions import *
from back_process import *
from time_integrators import *


if __name__ == '__main__':
    x_grid = np.arange(-100, 100, 0.1)
    dx = x_grid[1] - x_grid[0]
    A = 1.2
    mus = np.arange(-6, -0.01, 0.1)
    w = np.sqrt(2)
    X = np.arange(-5, 5, 0.05)
    F = []
    dF = []
    for j in range(len(mus)):
        Fj = []
        for i in range(len(X)):
            P = (1 + A * np.tanh(X[i]) ** 2)
            phi_01_i = np.sqrt(-2 * mus[j] / P) * (1 / np.cosh(np.sqrt(-mus[j]) * (x_grid - X[i])))
            Dphi_01_i = np.append(np.diff(phi_01_i) / (dx), 0)
            f_i = Dphi_01_i ** 2 \
                  + 0.5 * w ** 2 * x_grid ** 2 * phi_01_i ** 2 \
                  - 0.5 * (1 + A * np.tanh(x_grid) ** 2) * phi_01_i ** 4 - mus[j] * phi_01_i ** 2

            F_i = integrate.simpson(f_i, x_grid)
            Fj.append(F_i)
        Fj = np.array(Fj)
        dFj = np.diff(Fj) / (X[1] - X[0])
        F.append(Fj)
        dF.append(dFj)
    F = np.array(F)
    dF = np.array(dF)
    X, Y = np.meshgrid(X[:-1], mus)
    '''
    fig, (ax1) = plt.subplots(nrows=1, ncols=1)
    pc_01a = ax1.pcolor(X[1:], mus, np.abs(dF), cmap=parula_map)
    ax1.set_ylabel("$d$", fontsize=15)
    ax1.tick_params(axis="x", direction="in", labeltop=True, labelbottom=False, top=True, bottom=True)
    ax1.tick_params(axis="y", direction="in", left=True, right=True)
    ax1.grid(alpha=0.2, color="k")
    plt.show()
    plt.close()


    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    surf = ax.plot_surface(X, Y, np.abs(dF), cmap=parula_map, linewidth=0, antialiased=False)
    ax.set_xlabel("$b$")
    ax.set_ylabel("$A$")
    plt.show()
    plt.close()
    '''
    plt.contour(Y, X, dF, [0])
    plt.show()