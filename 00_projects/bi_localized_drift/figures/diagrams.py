import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    alpha = 1
    nu = 0.2
    sigma = 16
    d = 70
    k = np.sqrt(nu / alpha)
    x_grid = np.arange(-100, 100, 0.1)
    A = np.exp(-(x_grid + d / 2) ** 2 / (2 * sigma ** 2))
    B = np.exp(-(x_grid - d / 2) ** 2 / (2 * sigma ** 2))
    f = (A + B) * np.cos(k * x_grid)

    # Rango de parámetros
    param_min = 60
    param_max = 100

    # Tamaños para presentación
    label_fontsize = 28
    tick_fontsize = 20
    title_fontsize = 18
    point_size = 60

    fig, ax01 = plt.subplots(1, 1, figsize=(8, 3))

    # --- Gráfico superior izquierdo ---
    ax01.plot(x_grid, f, color="k", lw=3)
    ax01.plot(x_grid, A, color="purple", lw=3, ls="--")
    ax01.plot(x_grid, B, color="green", lw=3, ls="--")
    ax01.axis('off')
    plt.savefig("diagram01.png", dpi=300)