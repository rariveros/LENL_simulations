from functions import *
from back_process import *
from time_integrators import *
#from directories_lyap import *

if __name__ == '__main__':
    x_grid = np.arange(-250, 250)
    gamma_r = 1
    gamma_l = - 0 * gamma_r
    phi = 0 * np.pi
    dist = 1 * 33
    sigma = 7.9
    y_grid = gamma_l * np.exp(- (x_grid - dist / 2) ** 2 / (2 * sigma ** 2)) + gamma_r * np.cos(phi) * np.exp(- (x_grid + dist / 2) ** 2 / (2 * sigma ** 2))

    fig, ax = plt.subplots()
    ax.plot(x_grid, y_grid, "k")
    ax.set_aspect(20)
    plt.xlim([-150, 150])
    plt.ylim([-1.1, 1.1])
    plt.grid(alpha=0.3)
    plt.tick_params(left=False, right=False, labelleft=False, labelbottom=False, bottom=False)
    plt.savefig("forcing_03.png", dpi=200)