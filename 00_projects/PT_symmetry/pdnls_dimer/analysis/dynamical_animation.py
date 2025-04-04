import matplotlib.pyplot as plt
from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    directory = "D:/mnustes_science/simulation_data/FD/pdnlS_dimer/nu=0.0000/mu=0.1000/k=0.0500/gamma=0.1200"
    U1 = np.loadtxt(directory + '/U1.txt', delimiter=',', dtype="complex")
    V1 = np.loadtxt(directory + '/V1.txt', delimiter=',', dtype="complex")
    t = np.loadtxt(directory + '/T.txt', delimiter=',')

    # Crear figura y ejes
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(9, 4), facecolor='black')
    ax1.set_xlim([-0.4, 0.4])
    ax1.set_ylim([-0.4, 0.4])
    ax1.tick_params(colors='black', labelsize=18)  # color de las marcas y números
    ax1.set_xlabel("$Re\ \psi_1$", fontsize=20)
    ax1.set_ylabel("$Im\ \psi_1$", fontsize=20)

    ax2.set_xlim([-0.4, 0.4])
    ax2.set_ylim([-0.4, 0.4])
    ax2.tick_params(colors='black', labelsize=18)  # color de las marcas y números
    ax2.set_xlabel("$Re\ \psi_2$", fontsize=20)
    ax2.set_ylabel("$Im\ \psi_2$", fontsize=20)

    each = 5
    until = int(len(t)*(140 /1000))
    R1 = np.real(U1)[:until:each]
    I1 = np.imag(U1)[:until:each]
    R2 = np.real(V1)[:until:each]
    I2 = np.imag(V1)[:until:each]
    t = t[:until-1:each]
    print(len(t))
    print(len(R1))
    # Punto animado
    scat_1 = ax1.scatter(R1[0], I2[0], s=10, c="k", zorder=10)
    scat_2 = ax2.scatter(R2[0], I1[0], s=10, c="k", zorder=10)
    plt.tight_layout()

    # Función de animación
    def animate(i):
        scat_1.set_offsets((R1[i], I2[i]))
        scat_2.set_offsets((R2[i], I1[i]))
        return scat_1, scat_2

    # Animar
    ani = animation.FuncAnimation(fig, animate, frames=len(t), interval=50, blit=True, repeat=True)

    # Guardar como GIF
    ani.save("rabi.gif", dpi=250, writer=PillowWriter(fps=30))