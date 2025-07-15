from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':

    # Parámetros físicos
    L = 3.0  # tamaño del dominio [-L, L]
    Nx = 901
    Ny = 901
    dx = 2 * L / (Nx - 1)
    dy = dx  # malla cuadrada

    x = np.linspace(-L, L, Nx)
    y = np.linspace(-L, L, Ny)
    X, Y = np.meshgrid(x, y)

    # Velocidad en función de x (con salto en x=0)
    c1 = 2.0
    c2 = 4.0
    c = np.where(X < 0, c1, c2)

    # Parámetros del pulso
    x0, y0 = -1.5, -1.5  # posición inicial del centro del pulso
    sigma = 0.04  # ancho del pulso
    theta = np.radians(45)  # ángulo de incidencia (grados a radianes)
    lam = 0.12  # longitud de onda
    k = 2 * np.pi / lam

    # Componentes del vector de onda
    kx = k * np.cos(theta)
    ky = k * np.sin(theta)

    # Definir los vectores director y transversal
    n = np.array([np.cos(theta), np.sin(theta)])
    t = np.array([-np.sin(theta), np.cos(theta)])

    # Calcular las coordenadas en el sistema rotado
    X_shift = X - x0
    Y_shift = Y - y0
    s = n[0] * X_shift + n[1] * Y_shift
    tau = t[0] * X_shift + t[1] * Y_shift

    # Parámetros del haz
    sigma_s = 0.07  # longitud en la dirección de propagación
    sigma_t = 0.2  # ancho transversal (grande para evitar emisión lateral)
    s0 = 0.8

    # Condición inicial tipo haz gaussiano
    u0 = np.exp(-((s - s0) ** 2) / (2 * sigma_s ** 2)) * np.exp(-(tau ** 2) / (2 * sigma_t ** 2)) * np.cos(k * (s - s0))
    u = u0.copy()
    u_prev = u0.copy()

    # Tiempo
    c_max = max(c1, c2)
    dt = (0.4 * dx / c_max)  # condición CFL
    Nt = 2000

    # Generar capa absorbente (sponge layer)
    damping = np.ones_like(u)
    width = 1  # número de puntos de la capa absorbente
    strength = 0.01  # intensidad de la absorción

    for i in range(Nx):
        for j in range(Ny):
            ix = min(i, Nx - 1 - i)
            jy = min(j, Ny - 1 - j)
            dist = min(ix, jy)
            if dist < width:
                damping[i, j] = np.exp(-strength * (width - dist) ** 2)

    # Preparar animación
    fig, ax = plt.subplots()
    im = ax.imshow(u, extent=(-L, L, -L, L), cmap=parula_map, vmin=-0.25, vmax=0.25, animated=True)
    plt.colorbar(im)
    ax.axvline(0, color='black', linewidth=1)
    ax.set_xlim(-1.0, 1.0)
    ax.set_ylim(-1.0, 1.0)
    ax.set_xlabel('$x\ (\\textrm{mm})$')
    ax.set_ylabel('$y\ (\\textrm{mm})$')

    frames = []

    for n in range(Nt):
        laplacian = (
                            np.roll(u, 1, axis=0) + np.roll(u, -1, axis=0)
                            + np.roll(u, 1, axis=1) + np.roll(u, -1, axis=1)
                            - 4 * u
                    ) / dx ** 2

        u_new = 2 * u - u_prev + (dt ** 2) * c ** 2 * laplacian

        # Aplicar el amortiguamiento
        u_new *= damping

        u_prev = u
        u = u_new

        if n % 5 == 0:
            frame = ax.imshow(u, extent=(-L, L, -L, L), cmap=parula_map, vmin=-0.25, vmax=0.25, animated=True)
            ax.axvline(0, color='black', linewidth=1)
            frames.append([frame])

    # Crear la animación
    ani = animation.ArtistAnimation(fig, frames, interval=50, blit=True)
    ani.save("refraccion_onda_angular.gif", writer="pillow", dpi=300)