import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    # Datos de ejemplo
    x = np.linspace(0, 1, 500)
    k = 2 * np.pi / 1    # número de onda
    w = 2 * np.pi / 1    # frecuencia angular

    # Elongación de la onda
    y = np.sin(k * x)

    # Crear figura
    fig, ax = plt.subplots(figsize=(6, 3))

    # Graficar la línea
    ax.plot(x, y, color='#393534', lw=8)
    ax.set_ylim(-1.2, 1.2)

    ax.axis('off')
    ax.set_facecolor('white')  # fondo blanco limpio
    fig.patch.set_facecolor('white')

    plt.tight_layout()
    plt.savefig("ciclo.png", dpi=100)