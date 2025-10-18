from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    mu = 0.075
    gamma = 0.181
    alpha = 6.524
    sigma = 15
    x_grid = np.arange(-20, 20, 0.1)
    dx = x_grid[1] - x_grid[0]
    nus = np.arange(-10, 0.0, 0.1)
    G = []
    for i in range(len(nus)):
        nu = nus[i]
        delta = - nu + np.sqrt(gamma ** 2 - mu ** 2)
        k = np.sqrt(np.pi / 2) #np.arccosh(2)/(np.sqrt(2 * np.log(2)))#(np.sqrt(np.pi / 2)) #1.68 #(np.sqrt(np.pi / 2))
        print(k)
        sig_eff = k * np.sqrt(alpha / delta)  # np.arccosh(2)/(np.sqrt(2 * np.log(2))) * np.sqrt(alpha / delta)
        gamma_eff = gamma / np.sqrt(1 + sig_eff ** 2 / (sigma ** 2))
        #delta = - nu + np.sqrt(gamma_eff ** 2 - mu ** 2)
        lambd = np.sqrt(2 / np.pi) * sigma ** 2 * k / (sigma ** 2 + sig_eff ** 2) #np.sqrt(2 * delta) * sigma ** 2 / (sigma ** 2 + np.pi * alpha / (2 * delta))
        G_i = - mu * x_grid + mu * lambd * x_grid * np.exp(- x_grid ** 2 /  (2 * (sigma ** 2 + sig_eff ** 2)))
        G.append(G_i)
    G = np.array(G)

    pcm = plt.pcolormesh(x_grid, nus, G, cmap="seismic", shading='auto') # , vmin=-0.1, vmax=0.1
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label(r'$F(\xi)$', rotation=0, size=20, labelpad=-20, y=1.1)
    cbar.ax.tick_params(labelsize=15)
    cont = plt.contour(x_grid, nus, G, levels=[0], colors='k', linewidths=3)
    plt.gca().invert_yaxis()
    plt.xlim([x_grid[0], x_grid[-1]])
    plt.xlabel('$x$', size=25)
    plt.ylabel(r'$\nu$', size=25)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.grid(linestyle='--', alpha=0.2, color='k')
    plt.savefig('naive_bifurcation.svg', rasterized=True, dpi=300)
    plt.close()