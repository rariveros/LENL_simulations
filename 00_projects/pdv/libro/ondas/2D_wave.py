import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.mplot3d import Axes3D

if __name__ == '__main__':
    # Colormap azul-blanco-amarillo
    colors = [(0, '#0000ff'), (0.5, 'white'), (1, '#ffc621')]
    custom_cmap = LinearSegmentedColormap.from_list('custom', colors)

    # Dominio espacial
    r_max = 10
    x = np.linspace(-r_max, r_max, 500)
    y = np.linspace(-r_max, r_max, 500)
    X, Y = np.meshgrid(x, y)
    R = np.sqrt(X**2 + Y**2)

    # Parámetros de la onda radial
    k = 2 * np.pi / 2  # longitud de onda = 2
    w = 2 * np.pi / 1  # frecuencia = 1
    t = 0  # instante fijo

    # Campo de onda radial con decaimiento suave
    U = np.cos(k * R - w * t) * np.exp(-0.2 * R)

    # Ocultar todo lo que esté fuera de la región circular R > r_max
    U[R > r_max] = np.nan

    # Crear figura y eje 3D
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection='3d')

    # Dibujar superficie
    surf = ax.plot_surface(
        X, Y, U,
        cmap=custom_cmap,
        edgecolor='k',
        linewidth=0.0,
        antialiased=True,
        rstride=1, cstride=1,
        alpha=1.0,
        vmin=-0.9, vmax=0.9 # alpha original suave
    )

    # Estética limpia
    ax.set_zlim(-5, 5)
    ax.set_xlim(-r_max, r_max)
    ax.set_ylim(-r_max, r_max)
    ax.set_axis_off()
    ax.grid(False)
    ax.set_facecolor('white')
    fig.patch.set_facecolor('white')
    ax.view_init(elev=35, azim=135)

    # Mostrar
    plt.savefig("2D_wave.png", dpi=300, bbox_inches='tight')