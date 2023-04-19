from directories import *
from back_process import *


if __name__ == "__main__":

    disco = 'F'
    initial_dir_data = str(disco) + ':/mnustes_science/experimental_data'
    dirs = []
    N_dirs = 2

    x_ticks = [-20, -15, -10, -5, 0, 5, 10, 15, 20]
    t_i = 0
    t_f = -1

    for i in range(N_dirs):
        root = tk.Tk()
        root.withdraw()
        directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
        dirs.append(directory)
        initial_dir_data = os.path.dirname(directory)

    fig = plt.figure()
    fig.set_figheight(5)
    fig.set_figwidth(9)
    ax1 = plt.subplot2grid(shape=(7, 17), loc=(1, 0), colspan=8, rowspan=5)
    ax2 = plt.subplot2grid(shape=(7, 17), loc=(1, 8), colspan=8, rowspan=5)
    cbar_ax = plt.subplot2grid(shape=(7, 18), loc=(1, 17), colspan=1, rowspan=5)
    ax = [ax1, ax2]

    for i in range(N_dirs):
        X = np.loadtxt(dirs[i] + '/X.txt', delimiter=',')
        T = np.loadtxt(dirs[i] + '/T.txt', delimiter=',')
        Z = np.loadtxt(dirs[i] + '/field_real.txt', delimiter=',')

        norm = TwoSlopeNorm(vmin=np.amin(Z), vcenter=0, vmax=np.amax(Z))
        pcm = ax[i].pcolormesh(X, T[t_i:t_f], Z[t_i:t_f, :], norm=norm, cmap='seismic', shading='auto')

        ax[i].set_xlim([X[0], X[-1]])
        ax[i].set_xlabel('$x\ (\\textrm{mm})$', size='15')
        if i == 0:
            ax[i].set_ylabel('$t\ (\\textrm{s})$', size='15')
        else:
            plt.setp(ax[i].get_yticklabels(), visible=False)
            ax[i].tick_params(axis='y', which='both', length=0)
        ax[i].tick_params(axis='x', which='major', labelsize=8)
        ax[i].set_xticks(x_ticks)
        #ax[i].grid(linestyle='--', alpha=0.5)

    fig.subplots_adjust(left=0.1,
                        bottom=0.1,
                        right=0.9,
                        top=0.9,
                        wspace=1.2,
                        hspace=1.2)
    cbar = fig.colorbar(pcm, cax=cbar_ax)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label('$A_{R}(x, t)$', rotation=0, size=30, labelpad=-25, y=1.1)
    #plt.show()
    plt.savefig('E:/mnustes_science/simulation_data/tesis_rafael/pattern_formation.png', dpi=300)
    plt.close()