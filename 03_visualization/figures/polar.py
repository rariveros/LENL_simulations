from directories import *
from back_process import *

if __name__ == '__main__':
    r = np.arange(0, 2, 0.01)
    t = np.arange(0, 2 * np.pi, 0.1)
    Nt = len(t)
    Nr = len(r)

    R, T = np.meshgrid(r, t)

    theta = T
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.plot(theta, r, lw = 3)
    ax.set_rmax(2)
    ax.grid(False)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    plt.show()