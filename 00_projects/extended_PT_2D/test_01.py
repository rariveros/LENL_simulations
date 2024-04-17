from back_process import *
from functions import *
from time_integrators import *


from functions import *
from back_process import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    N = 8
    Lx = 100
    Ly = 100
    dx = 0.75
    dy = 0.75
    sigma = 6

    [xmin, xmax, dx] = [- Lx / 2,  Lx / 2, dx]
    [ymin, ymax, dy] = [- Ly / 2, Ly / 2, dy]
    x_grid = np.arange(xmin, xmax, dx)
    y_grid = np.arange(ymin, ymax, dy)
    X, Y = np.meshgrid(x_grid, y_grid)
    positions = []
    Zs = []
    R = 20
    for i in range(N):
        x_i = R * np.cos(2 * np.pi * i / N)
        y_i = R * np.sin(2 * np.pi * i / N)
        Z_i = ((-1) ** i) * np.exp(-((X - x_i) ** 2 + (Y - y_i) ** 2) / (2 * sigma ** 2))
        positions.append([x_i, y_i])
        Zs.append(Z_i)
    Zs = np.array(Zs)
    Z = Zs[0]
    for i in range(1, len(Zs)):
        Z = Z + Zs[i]
    positions = np.array(positions)
    pcm = plt.pcolormesh(x_grid, y_grid, Z, cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$\gamma(\\vec{x})$', rotation=0, size=20, labelpad=-27, y=1.15)
    plt.xlabel('$y$', size='20')
    plt.ylabel('$x$', size='20')
    plt.grid(linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.show()
    plt.close()

    plt.scatter(positions[:, 0], positions[:, 1], zorder=10)
    plt.scatter(np.mean(positions[:, 0]), np.mean(positions[:, 1]), zorder=10)
    plt.xlim([-1.1, 1.1])
    plt.ylim([-1.1, 1.1])
    plt.grid(alpha=0.35)
    plt.show()

    ### FUNCIONA ###