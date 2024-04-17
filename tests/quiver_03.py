from functions import *
from back_process import *
from time_integrators import *


if __name__ == '__main__':
    # Creating arrow
    x_min = -0.5
    x_max = - x_min
    y_min = -0.5
    y_max = - y_min
    x = np.arange(x_min, x_max, x_max/1000)
    y = np.arange(y_min, y_max, y_max/1000)
    z = [0]
    x_grid = np.arange(x_min, x_max, x_max/1000)

    X, Y = np.meshgrid(x, y)
    #X, Y, Z = np.meshgrid(x, y, z)
    nu = 0.25
    mu = 0.1
    gamma = 0.2
    beta = 1
    u = (beta * (X ** 2 + Y ** 2) - nu) * Y + (gamma - mu) * X
    v = - (beta * (X ** 2 + Y ** 2) - nu) * X - (gamma + mu) * Y
    #w = X + Y - 2 * Z + X ** 2 * Y

    # creating plot
    #ax = plt.figure().add_subplot(projection='3d')
    #ax.quiver(X, Y, Z, u, v, w, length=0.05, normalize=True, arrow_length_ratio=0.3)

    #fig, axs = plt.subplots(1, 1)
    #ax1 = axs
    #q_0 = ax1.quiver(X, Y, u / np.sqrt(u ** 2 + v ** 2), v / np.sqrt(u ** 2 + v ** 2), scale=25, alpha=0.6)
    #q_0 = ax1.quiver(X, Y, u, v, scale=1, alpha=0.6)
    #ax1.scatter([0], [0], c="r")
    #ax1.axis([x_min + 0.1, x_max - 0.1, y_min + 0.1, y_max - 0.1])
    #ax1.set_aspect("equal")
    #plt.grid(alpha=0.3)
    #ax1.set_title('$\mu = 0.4$')
    #ax1.set_xlabel('$x$', fontsize=12)
    #ax1.set_ylabel('$y$', fontsize=12)
    #plt.show()
    #plt.close()
    #plt.tight_layout()
    #plt.savefig("hopf.png", dpi=300)

    #fig, (ax1) = plt.subplots(nrows=1, ncols=1)

    #pc_01 = ax1.pcolor(x, y, np.log(u ** 2 + v ** 2), cmap=parula_map)
    #ax1.grid(alpha=0.2, color="k")
    #ax1.tick_params(axis="x", direction="in", top=True, bottom=True)
    #ax1.tick_params(axis="y", direction="in", left=True, right=True)
    #plt.show()
    #plt.close()

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    surf = ax.plot_surface(X, Y, np.log(u ** 2 + v ** 2), cmap=parula_map, linewidth=0, antialiased=False)
    plt.show()