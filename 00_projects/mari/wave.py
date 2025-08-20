from back_process import *
from matplotlib.patches import Circle

if __name__ == '__main__':
    XL = np.arange(-6, 0.02, 0.02)
    XR = np.arange(0, 6, 0.02)
    t = np.arange(0, 2 * np.pi, 0.05)
    F1 = np.exp(np.pi * 1j * XL)
    F2 = np.exp(np.pi * 1j * XR)  # Crear figura y ejes
    X0L = -4
    X0R = 4
    I0L = np.argmin(np.abs(XL - X0L))
    I0R = np.argmin(np.abs(XR - X0R))

    fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(5, 2))#, facecolor='black')
    ax1.set_xlim([-3, 3])
    ax1.set_ylim([-1.2, 1.2])
    #ax1.set_facecolor('black')
    ax1.tick_params(colors='k', labelsize=18)  # color de las marcas y números
    #ax1.spines['bottom'].set_color('white')
    #ax1.spines['top'].set_color('white')
    #ax1.spines['right'].set_color('white')
    #ax1.spines['left'].set_color('white')
    #ax1.yaxis.label.set_color('white')
    #ax1.xaxis.label.set_color('white')
    ax1.set_xlabel("$x\ \\textrm{(mm)}$", fontsize=20)
    #ax1.set_ylabel("$\Delta M$", fontsize=20)

    # Punto animado
    line_1, = ax1.plot(XL, np.real(F1 * np.exp(-1j * t[0])), color="k", zorder=10, lw=8)
    line_2, = ax1.plot(XR, np.real(F2 * np.exp(-1j * t[0])), color="k", zorder=10, lw=8)
    #scat_1 = ax1.scatter(X0L, np.real(F1[I0L] * np.exp(-1j * t[0])), s=80, c="k", zorder=10)
    #scat_2 = ax1.scatter(X0L, np.real(F2[I0R] * np.exp(-1j * t[0])), s=80, c="k", zorder=10)
    #scat_3 = ax1.scatter(X0L, np.real(F2[0] * np.exp(-1j * t[0])), s=80, c="k", zorder=10)
    #ax1.hlines(0, -6, 6, colors="k", linestyles="dashed", lw=1)
    plt.tight_layout()

    # Función de animación
    def animate(i):
        line_1.set_data(XL, np.real(F1 * np.exp(-1j * t[i + 1])))
        line_2.set_data(XR, np.real(F2 * np.exp(-1j * t[i + 1])))
        #scat_1.set_offsets((X0L, np.real(F1[I0L] * np.exp(-1j * t[i + 1]))))
        #scat_2.set_offsets((X0R, np.real(F2[I0R] * np.exp(-1j * t[i + 1]))))
        #scat_3.set_offsets((0, np.real(F2[0] * np.exp(-1j * t[i + 1]))))
        return line_1, line_2,# scat_1, scat_2, scat_3,

    # Animar
    ani = animation.FuncAnimation(fig, animate, frames=len(t) - 1, interval=50, blit=True, repeat=True)

    # Guardar como GIF
    ani.save("wave.gif", dpi=250, writer=PillowWriter(fps=30))