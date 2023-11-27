from back_process import *

if __name__ == '__main__':
    Lx = Ly = 20
    x = np.arange(-Lx / 2, Lx / 2, 0.05)
    y = np.arange(-10, 10, 0.05)
    t = np.arange(0, 10, 0.1)
    Nx = len(x)
    Ny = len(y)
    Nt = len(t)

    k_x = 1
    k_y = 1
    w = 1

    X, Y, T = np.meshgrid(x, y, t)
    XYT = np.cos(k_x * X) * np.cos(k_y * Y) * np.sin(w * T)

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.set(xlim=(-Lx / 2, Lx / 2), ylim=(-Ly / 2, Ly / 2))
    cax = ax.pcolormesh(x, y, XYT[:, :, 0], vmin=-1, vmax=1, cmap='RdBu', shading='auto')
    cbar = fig.colorbar(cax)
    cbar.set_label('$u(x, y)$', rotation=0, size=20, labelpad=-27, y=1.13)


    def animate(i):
        cax.set_array(XYT[:, :, i].flatten())


    anim = FuncAnimation(fig, animate, interval=1, frames=Nt - 1)

    FFwriter = animation.FFMpegWriter()
    anim.save('name_test.mp4', writer=FFwriter, dpi=300)
    plt.close()

    pcm = plt.pcolormesh(x, y, XYT[:, :, 0], vmin=-1, vmax=1, cmap='RdBu', shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$u(x, y)$', rotation=0, size=20, labelpad=-27, y=1.15)
    plt.xlabel('$y$', size='20')
    plt.ylabel('$x$', size='20')
    plt.grid(linestyle='--', alpha=0.5)
    #plt.savefig("TLI.gif", writer = 'imagemagick', fps = 60)
    #plt.savefig('final_profile.png', dpi=300)
    plt.close()
