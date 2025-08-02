import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

if __name__ == '__main__':
    # Dominio espacial y temporal
    x1 = np.linspace(0, 4.9, 250)
    x2 = np.linspace(5.05, 10, 250)
    x = np.linspace(0, 10, 500)
    t = np.linspace(0, 5, 500)
    X1, T = np.meshgrid(x1, t)
    X2, T = np.meshgrid(x2, t)

    # Parámetros de la onda
    k = 2 * np.pi / 3    # número de onda
    w = 2 * np.pi / 1    # frecuencia angular

    # Elongación de la onda
    U1 = np.sin(k * X1 - w * T)
    U2 = np.sin(k * X2 - w * T)

    # Gráfico
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection='3d')

    # Superficie suave y continua
    surf1 = ax.plot_surface(
        X1, T, U1,
        cmap='cividis',     # sin bordes
        edgecolor='k',  # agrega bordes
        linewidth=0.01,  # sin bordes
        antialiased=True,
        rstride=1, cstride=1
    )
    surf2 = ax.plot_surface(
        X2, T, U2,
        cmap='cividis',        # paleta continua amigable
        edgecolor='k',  # agrega bordes
        linewidth=0.01,  # sin bordes
        antialiased=True,
        rstride=1, cstride=1
    )
    ax.plot([5] * len(t), t, np.sin(k * 5 - w * t), color='red', lw=5, zorder=3)
    ax.plot(x, [5.02] * len(x), np.sin(k * x - w * 0), color='blue', lw=5, zorder=6)
    # Estética limpia
    ax.set_zlim(-10, 10)
    ax.set_axis_off()
    ax.grid(False)
    ax.set_facecolor('white')
    fig.patch.set_facecolor('white')

    # Ángulo de vista sugerido
    ax.view_init(elev=25, azim=120)

    plt.tight_layout()
    plt.savefig("espaciotemporal.png", dpi=300)
