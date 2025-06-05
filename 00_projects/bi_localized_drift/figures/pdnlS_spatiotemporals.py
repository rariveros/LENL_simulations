from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    frequencies = []
    disc = "aD:"
    initial_dir_data = str(disc) + 'mnustes_science/simulation_data/FD'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    x_grid = np.loadtxt(directory + '/X.txt', delimiter=',')
    t_light = np.loadtxt(directory + '/T.txt', delimiter=',')
    ZR = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    ZI = np.loadtxt(directory+ '/field_img.txt', delimiter=',')
    params = np.loadtxt(directory+ '/parameters.txt', delimiter=',')
    [alpha, beta, gamma_0, dist, sigma, mu, nu] = params

    Z = ZR + 1j * ZI
    modulo_light_1 = np.abs(Z)
    arg_light_1 = np.angle(Z)
    arg_light_1 = (2 * np.pi + arg_light_1) * (arg_light_1 < 0) + arg_light_1 * (arg_light_1 > 0)

    phi = np.pi
    gamma_complex = gamma_0 * (np.exp(- (x_grid - dist / 2) ** 2 / (2 * sigma ** 2)) + np.exp(1j * phi) * np.exp( - (x_grid + dist / 2) ** 2 / (2 * sigma ** 2)))
    gamma_real = np.real(gamma_complex)

    beta = 0.004811649356064012

    figsize = (5, 4)
    dpi = 300
    fontsize_labels = 24
    fontsize_ticks = 20
    x_min = x_grid[0]
    x_max = x_grid[-1]
    t_0 = 2000
    t_light = t_light - t_0
    t_min = 0
    t_max = 3001 #t_light[-1]
    grid_style = '--'
    grid_alpha = 0.0

    plt.figure(figsize=figsize)
    pcm = plt.pcolormesh(x_grid, t_light, modulo_light_1 / np.sqrt(beta), cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$|\psi|$', rotation=0, size=fontsize_labels, labelpad=-25, y=1.15)
    cbar.ax.tick_params(labelsize=fontsize_ticks)
    plt.xlabel('$x\ \\textrm{(mm)}$', fontsize=fontsize_labels)
    plt.ylabel('$t$', fontsize=fontsize_labels)
    plt.tick_params(labelsize=fontsize_ticks)
    plt.xlim([x_min, x_max])
    plt.ylim([t_min, t_max])
    plt.grid(linestyle=grid_style, alpha=grid_alpha)
    plt.tight_layout()
    plt.savefig(directory + '/module_pretty.png', dpi=dpi)
    plt.close()

    plt.figure(figsize=figsize)
    pcm = plt.pcolormesh(x_grid, t_light, arg_light_1, cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$\\textrm{arg}(\psi)$', rotation=0, size=fontsize_labels, labelpad=-25, y=1.15)
    cbar.ax.tick_params(labelsize=fontsize_ticks)
    plt.xlabel('$\ \\textrm{(mm)}$', fontsize=fontsize_labels)
    plt.ylabel('$t$', fontsize=fontsize_labels)
    plt.tick_params(labelsize=fontsize_ticks)
    plt.xlim([x_min, x_max])
    plt.ylim([t_min, t_max])
    plt.grid(linestyle=grid_style, alpha=grid_alpha)
    plt.tight_layout()
    plt.savefig(directory + '/arg_pretty.png', dpi=dpi)
    plt.close()

    fig, ax01 = plt.subplots(1, 1, figsize=(8, 2))
    #ax01.plot(x_grid, np.amax(modulo_light_1 / np.sqrt(beta)) * (gamma_real / gamma_0), color="k", label="$\gamma(x)$", lw=3)
    ax01.plot(x_grid, np.real(Z[0, :]/ np.sqrt(beta)), color="b", label="$\psi_R$", lw=5)
    ax01.plot(x_grid, np.imag(Z[0, :]/ np.sqrt(beta)), color="r", label="$\psi_I$", lw=5)
    ax01.set_xlim([x_grid[0], x_grid[-1]])
    ax01.tick_params(labelsize=fontsize_ticks, labelleft=False, labelright=False, left=False, right=False)
    #ax01.legend(loc="upper right", fontsize=20)
    plt.grid(linestyle='--', alpha=0.5)
    plt.savefig(directory + '/pretty_init_profiles.png', dpi=300, bbox_inches='tight')
    plt.close()