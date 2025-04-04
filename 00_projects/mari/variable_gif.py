from back_process import *
from matplotlib.patches import Circle

if __name__ == '__main__':
    T = 4
    t = np.arange(0, 2 * T + 0.025, 0.025)
    fact = 2 * np.pi
    # Función
    a1, a2, a3 = 5.2, 0.58, 0.37
    b1, b2, b3 = 3.4, -0.1, -0.44
    F = (a1 * np.cos(fact * t / T) + a2 * np.cos(2 * fact * t / T) + a3 * np.cos(3 * fact * t / T) + \
        b1 * np.sin(fact * t / T) + b2 * np.sin(2 * fact * t / T) + b3 * np.sin(3 * fact * t / T)) / 7

    width, height = 1000, 800

    # Número de estrellas
    num_stars = 300

    # Coordenadas aleatorias de las estrellas
    x = np.random.rand(num_stars) * width
    y = np.random.rand(num_stars) * height

    # Tamaños y brillos aleatorios
    sizes = np.random.rand(num_stars) * 3 + 0.5
    alphas = np.random.rand(num_stars) * 0.25 + 0.1

    # Colores RGBA (blanco con alfa variable)
    colors = np.ones((num_stars, 4))  # RGBA
    colors[:, 0:3] = 1  # blanco
    colors[:, 3] = alphas  # transparencia variable

    # Crear figura y ejes
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(9, 4), facecolor='black')
    ax1.set_xlim([0, 2 * T])
    ax1.set_ylim([-1, 1])
    ax1.set_facecolor('black')
    ax1.tick_params(colors='white', labelsize=18)  # color de las marcas y números
    ax1.spines['bottom'].set_color('white')
    ax1.spines['top'].set_color('white')
    ax1.spines['right'].set_color('white')
    ax1.spines['left'].set_color('white')
    ax1.yaxis.label.set_color('white')
    ax1.xaxis.label.set_color('white')
    ax1.set_xlabel("$\\textrm{Days}$", fontsize=20)
    ax1.set_ylabel("$\Delta M$", fontsize=20)

    # Punto animado
    scat_1 = ax1.scatter(t[0], F[0], s=300, c="r", zorder=10)
    ax1.plot(t, F, lw=2, color="white")

    circle = Circle((width / 2, height / 2), radius=200, fill=True, edgecolor='white', facecolor="white", linewidth=2)
    ax2.add_patch(circle)
    ax2.set_facecolor('black')
    ax2.scatter(x, y, s=sizes, color=colors, marker='*')
    ax2.set_xlim(0, width)
    ax2.set_ylim(0, height)
    ax2.tick_params(colors='white', labelsize=20)
    ax2.set_xticks([])
    ax2.set_yticks([])
    # color de las marcas y números
    ax2.spines['bottom'].set_color('white')
    ax2.spines['top'].set_color('white')
    ax2.spines['right'].set_color('white')
    ax2.spines['left'].set_color('white')
    ax2.yaxis.label.set_color('white')
    ax2.xaxis.label.set_color('white')

    plt.tight_layout()

    each_time = 2
    t_grid = t[::each_time]
    G = F[::each_time]
    # Función de animación
    def animate(i):
        scat_1.set_offsets((t_grid[i], G[i]))
        base_radius = 200
        circle.set_radius(base_radius + G[i] * 50)
        return scat_1,

    # Animar
    ani = animation.FuncAnimation(fig, animate, frames=len(t_grid), interval=50, blit=True, repeat=True)

    # Guardar como GIF
    ani.save("unbroken_PT.gif", dpi=250, writer=PillowWriter(fps=30))