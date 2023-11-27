from functions import *
from back_process import *
from time_integrators import *


if __name__ == '__main__':
    x_grid = np.arange(0, 100, 0.1)
    y_grid = np.full_like(x_grid, 0)
    y_scatter = np.full_like(x_grid, np.nan)
    critical_value = 0
    for i in range(len(x_grid)):
        noise = 0.05 * (np.random.rand(1) - 0.5)
        if x_grid[i] > critical_value:
            y_grid[i] = 0.05 * (x_grid[i] - critical_value) ** 0.5
            if x_grid[i] > 8:
                y_scatter[i] = 0.05 * (x_grid[i] - critical_value) ** 0.5 + noise
    plt.scatter(x_grid[::50], y_scatter[::50], c="k", zorder=10)
    plt.plot(x_grid, y_grid, c="r", linestyle="--")
    plt.grid(alpha=0.3)
    plt.xlim([0, 100])
    plt.xlabel("$\sigma_i$", fontsize=20)
    plt.ylabel("$\sigma_w$", fontsize=20)
    plt.show()