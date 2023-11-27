from back_process import *


if __name__ == "__main__":
    disco = 'C:/'
    initial_dir_data = str(disco) + ':/mnustes_science/simulation_data'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    Z = np.loadtxt(directory + '/field.dat')

    pcm = plt.pcolormesh((np.arange(len(Z[0, :])) - len(Z[0, :]) / 2) * 0.25, np.arange(len(Z[:, 0])) / 10, Z / np.sqrt(0.00481), cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm)
    cbar.set_label('$|A(x,t)|$', rotation=0, size=25, labelpad=-27, y=1.11)

    # put the major ticks at the middle of each cell

    plt.xlim([-150, 150])
    plt.xlabel('$x$', size='25')
    plt.xticks(fontsize=15)

    plt.ylim([0, 610])
    plt.ylabel('$t$', size='25')
    plt.yticks(fontsize=15)

    plt.grid(linestyle='--', alpha=0.2, color='k')
    cbar.ax.tick_params(labelsize=15)

    # labels
    plt.tight_layout()
    #plt.show()
    plt.savefig('field_visual.png', dpi=300)
    plt.close()

    print(Z)