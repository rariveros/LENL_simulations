import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

if __name__ == '__main__':
    # Parámetros
    R_max = 1.0
    N = 300
    n = 5  # Para tener 2.5 ciclos → onda termina en máximo en R = 1
    k = (2 * n + 1) * np.pi / 2 / R_max

    # Malla cartesiana
    x = np.linspace(-R_max, R_max, N)
    y = np.linspace(-R_max, R_max, N)
    z = np.linspace(-R_max, R_max, N)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    # Coordenadas esféricas
    R = np.sqrt(X**2 + Y**2 + Z**2)
    PHI = np.arctan2(Y, X) % (2 * np.pi)

    # Onda esférica con decaimiento suave
    U = np.sin(k * R) * np.exp(-0.3 * R)

    # Máscara para cortar un gajo en phi ∈ [0, π/2]
    mask = ~((PHI >= 0) & (PHI <= 0.7 * np.pi)) & (R <= R_max)

    # Datos filtrados
    Xf, Yf, Zf = X[mask], Y[mask], Z[mask]
    Uf = U[mask]

    # Colormap rojo-blanco-azul
    colors = [(0, '#0000ff'), (0.5, 'white'), (1, '#ffc621')]
    rwb_cmap = LinearSegmentedColormap.from_list('rwb', colors)

    # Crear figura
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    sc = ax.scatter(Xf, Yf, Zf, c=Uf, cmap=rwb_cmap, s=1, vmin=-1, vmax=1)

    # Vista limpia y proporciones iguales
    ax.set_xlim(-R_max, R_max)
    ax.set_ylim(-R_max, R_max)
    ax.set_zlim(-R_max, R_max)
    ax.set_box_aspect([1, 1, 1])  # relación de aspecto 1:1:1
    ax.view_init(elev=25, azim=45)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.grid(False)
    ax.set_axis_off()

    # Colorbar más ancho y sin ticks
    cb = plt.colorbar(sc, ax=ax, shrink=0.8, pad=0.05, alpha=0.4)
    cb.ax.tick_params(labelsize=0, length=0)
    cb.outline.set_visible(False)

    plt.tight_layout()
    plt.savefig("3D_wave.png", dpi=300, bbox_inches='tight')