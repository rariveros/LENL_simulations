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
    A = np.loadtxt(directory + '/ampd_field_real.txt', delimiter=',', dtype=complex)
    B = np.loadtxt(directory+ '/ampd_field_img.txt', delimiter=',', dtype=complex)
    params = np.loadtxt(directory+ '/parameters.txt', delimiter=',')
    modulo_light_1 = np.abs(A + B)
    arg_light_1 = np.angle(A)
    arg_light_1 = (2 * np.pi + arg_light_1) * (arg_light_1 < 0) + arg_light_1 * (arg_light_1 > 0)
    arg_light_2 = np.angle(B)
    arg_light_2 = (2 * np.pi + arg_light_1) * (arg_light_1 < 0) + arg_light_1 * (arg_light_1 > 0)

    beta = 0.004811649356064012

    figsize = (5, 4)
    dpi = 300
    fontsize_labels = 20
    fontsize_ticks = 14
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
    cbar.set_label('$|A + B|$', rotation=0, size=fontsize_labels, labelpad=-25, y=1.15)
    cbar.ax.tick_params(labelsize=fontsize_ticks)
    plt.xlabel('$x\ \\textrm{(mm)}$', fontsize=fontsize_labels)
    plt.ylabel('$t$', fontsize=fontsize_labels)
    plt.tick_params(labelsize=fontsize_ticks)
    plt.xlim([x_min, x_max])
    plt.ylim([t_min, t_max])
    plt.grid(linestyle=grid_style, alpha=grid_alpha)
    plt.tight_layout()
    plt.savefig(directory + '/AB_module_pretty.png', dpi=dpi)
    plt.close()

    plt.figure(figsize=figsize)
    pcm = plt.pcolormesh(x_grid, t_light, np.abs(A) / np.sqrt(beta), cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$|A|$', rotation=0, size=fontsize_labels, labelpad=-25, y=1.15)
    cbar.ax.tick_params(labelsize=fontsize_ticks)
    plt.xlabel('$x\ \\textrm{(mm)}$', fontsize=fontsize_labels)
    plt.ylabel('$t$', fontsize=fontsize_labels)
    plt.tick_params(labelsize=fontsize_ticks)
    plt.xlim([x_min, x_max])
    plt.ylim([t_min, t_max])
    plt.grid(linestyle=grid_style, alpha=grid_alpha)
    plt.tight_layout()
    plt.savefig(directory + '/A_module_pretty.png', dpi=dpi)
    plt.close()

    plt.figure(figsize=figsize)
    pcm = plt.pcolormesh(x_grid, t_light, np.abs(B) / np.sqrt(beta), cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$|B|$', rotation=0, size=fontsize_labels, labelpad=-25, y=1.15)
    cbar.ax.tick_params(labelsize=fontsize_ticks)
    plt.xlabel('$x\ \\textrm{(mm)}$', fontsize=fontsize_labels)
    plt.ylabel('$t$', fontsize=fontsize_labels)
    plt.tick_params(labelsize=fontsize_ticks)
    plt.xlim([x_min, x_max])
    plt.ylim([t_min, t_max])
    plt.grid(linestyle=grid_style, alpha=grid_alpha)
    plt.tight_layout()
    plt.savefig(directory + '/B_module_pretty.png', dpi=dpi)
    plt.close()

    plt.figure(figsize=figsize)
    pcm = plt.pcolormesh(x_grid, t_light, arg_light_1, cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$\\textrm{arg}(A)$', rotation=0, size=fontsize_labels, labelpad=-25, y=1.15)
    cbar.ax.tick_params(labelsize=fontsize_ticks)
    plt.xlabel('$\ \\textrm{(mm)}$', fontsize=fontsize_labels)
    plt.ylabel('$t$', fontsize=fontsize_labels)
    plt.tick_params(labelsize=fontsize_ticks)
    plt.xlim([x_min, x_max])
    plt.ylim([t_min, t_max])
    plt.grid(linestyle=grid_style, alpha=grid_alpha)
    plt.tight_layout()
    plt.savefig(directory + '/A_arg_pretty.png', dpi=dpi)
    plt.close()

    plt.figure(figsize=figsize)
    pcm = plt.pcolormesh(x_grid, t_light, arg_light_2, cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$\\textrm{arg}(B)$', rotation=0, size=fontsize_labels, labelpad=-25, y=1.15)
    cbar.ax.tick_params(labelsize=fontsize_ticks)
    plt.xlabel('$\ \\textrm{(mm)}$', fontsize=fontsize_labels)
    plt.ylabel('$t$', fontsize=fontsize_labels)
    plt.tick_params(labelsize=fontsize_ticks)
    plt.xlim([x_min, x_max])
    plt.ylim([t_min, t_max])
    plt.grid(linestyle=grid_style, alpha=grid_alpha)
    plt.tight_layout()
    plt.savefig(directory + '/B_arg_pretty.png', dpi=dpi)
    plt.close()