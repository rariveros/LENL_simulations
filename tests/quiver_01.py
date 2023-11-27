from functions import *
from back_process import *
from time_integrators import *


if __name__ == '__main__':
    # Creating arrow
    x_min = -5
    x_max = - x_min
    y_min = -5
    y_max = - y_min
    x = np.arange(x_min, x_max, x_max/20)
    y = np.arange(y_min, y_max, y_max/20)
    z = [0]
    x_grid = np.arange(x_min, x_max, x_max/1000)
    mu_0 = 1
    mu_1 = 0
    mu_2 = -1
    a = 0
    X, Y = np.meshgrid(x, y)
    #X, Y, Z = np.meshgrid(x, y, z)
    u = Y
    v_0 = - mu_0 * (X ** 2 - 1) * Y - X + a
    v_1 = - mu_1 * (X ** 2 - 1) * Y - X + a
    v_2 = - mu_2 * (X ** 2 - 1) * Y - X + a
    #w = X + Y - 2 * Z + X ** 2 * Y

    # creating plot
    #ax = plt.figure().add_subplot(projection='3d')
    #ax.quiver(X, Y, Z, u, v, w, length=0.05, normalize=True, arrow_length_ratio=0.3)

    fig, axs = plt.subplots(1, 3)
    ax1 = axs[0]
    ax2 = axs[1]
    ax3 = axs[2]
    q_0 = ax1.quiver(X, Y, u, v_0, scale=150, alpha=0.6)
    ax1.scatter([a], [0], c="r")
    q_1 = ax2.quiver(X, Y, u, v_1, scale=100, alpha=0.6)
    ax2.scatter([a], [0], c="r")
    q_2 = ax3.quiver(X, Y, u, v_2, scale=150, alpha=0.6)
    ax3.scatter([a], [0], c="r")

    ax1.axis([x_min + 1, x_max - 1, y_min + 1, y_max - 1])
    ax1.set_aspect("equal")
    ax2.axis([x_min + 1, x_max - 1, y_min + 1, y_max - 1])
    ax2.set_aspect("equal")
    ax3.axis([x_min + 1, x_max - 1, y_min + 1, y_max - 1])
    ax3.set_aspect("equal")
    plt.grid(alpha=0.3)

    ax1.set_title('$\mu = 0.4$')
    ax1.set_xlabel('$x$', fontsize=12)
    ax1.set_ylabel('$y$', fontsize=12)

    ax2.set_title('$\mu = 0.0$')
    ax2.set_xlabel('$x$', fontsize=12)

    ax3.set_title('$\mu = -0.4$')
    ax3.set_xlabel('$x$', fontsize=12)

    plt.tight_layout()
    plt.savefig("hopf.png", dpi=300)