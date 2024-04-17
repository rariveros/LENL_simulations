from back_process import *

if __name__ == '__main__':
    Lx = 75
    Ly = 250
    x = np.arange(-Lx, Lx, 0.5)
    y = np.arange(0, Ly, 0.1)
    t = np.arange(0, 10, 0.1)
    Nx = len(x)
    Ny = len(y)
    Nt = len(t)

    n_x = 1
    n_y = 3
    w = 0.1

    k_x = n_x * np.pi / Lx
    k_y = n_y * np.pi / Ly

    X, Y = np.meshgrid(x, y)
    f = np.exp(- (X - 5) ** 2 / (2 * 6 ** 2))
    g = f + np.flip(f)
    h = f - np.flip(f)
    phi = g * np.exp(-1j * w * Y) + h * np.exp(1j * w * Y)
    v_x = k_x * np.sin(k_x * X) * np.cos(k_y * Y)
    v_y = k_y * np.cos(k_x * X) * np.sin(k_y * Y)

    pcm = plt.pcolormesh(x, y, np.real(phi), vmin=-2, vmax=2, cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$\phi_R(x, t)$', rotation=0, size=20, labelpad=-27, y=1.15)
    plt.xlabel('$x$', size='20')
    plt.ylabel('$t$', size='20')
    plt.grid(linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig('phi_R.png', dpi=150)
    plt.close()

    pcm = plt.pcolormesh(x, y, np.imag(phi), vmin=-2, vmax=2, cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$\phi_I(x, t)$', rotation=0, size=20, labelpad=-27, y=1.15)
    plt.xlabel('$x$', size='20')
    plt.ylabel('$t$', size='20')
    plt.grid(linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig('phi_I.png', dpi=150)
    plt.close()

    pcm = plt.pcolormesh(x, y, np.abs(phi), vmin=0, vmax=2, cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$|\phi(x, t)|$', rotation=0, size=20, labelpad=-27, y=1.15)
    plt.xlabel('$x$', size='20')
    plt.ylabel('$t$', size='20')
    plt.grid(linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig('phi_module.png', dpi=150)
    plt.close()

    skip = (slice(None, None, 50), slice(None, None, 50))

    plt.quiver(X[skip], Y[skip], v_x[skip], v_y[skip])
    plt.title("$\\textrm{Velocity Field}\ v(x,y)$", size="25")
    plt.xlim([x[0], x[-1]])
    plt.ylim([y[0], y[-1]])
    plt.xlabel('$y$', size='20')
    plt.ylabel('$x$', size='20')
    plt.tight_layout()
    plt.savefig('v_field.png', dpi=150)
    plt.close()