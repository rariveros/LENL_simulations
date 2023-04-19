from directories import *
from back_process import *

if __name__ == '__main__':
    Lx = Ly = 20
    x = np.arange(0, Lx, 0.01)
    y = np.arange(0, Ly, 0.01)
    t = np.arange(0, 10, 0.1)
    Nx = len(x)
    Ny = len(y)
    Nt = len(t)

    n_x = 1
    n_y = 3

    k_x = n_x * np.pi / Lx
    k_y = n_y * np.pi / Ly

    X, Y= np.meshgrid(x, y)
    phi = np.cos(k_x * X) * np.cos(k_y * Y)
    v_x = k_x * np.sin(k_x * X) * np.cos(k_y * Y)
    v_y = k_y * np.cos(k_x * X) * np.sin(k_y * Y)

    pcm = plt.pcolormesh(x, y, phi, vmin=-1, vmax=1, cmap='RdBu', shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$\phi(x, y)$', rotation=0, size=20, labelpad=-27, y=1.15)
    plt.xlabel('$y$', size='20')
    plt.ylabel('$x$', size='20')
    plt.grid(linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig('phi.png', dpi=300)
    plt.close()

    skip = (slice(None, None, 50), slice(None, None, 50))

    plt.quiver(X[skip], Y[skip], v_x[skip], v_y[skip])
    plt.title("$\\textrm{Velocity Field}\ v(x,y)$", size="25")
    plt.xlim([x[0], x[-1]])
    plt.ylim([y[0], y[-1]])
    plt.xlabel('$y$', size='20')
    plt.ylabel('$x$', size='20')
    plt.tight_layout()
    plt.savefig('v_field.png', dpi=300)
    plt.close()