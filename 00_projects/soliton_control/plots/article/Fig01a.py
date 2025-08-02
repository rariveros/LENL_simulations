from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    working_directory = "D:/mnustes_science/experimental_data/soliton_control"
    shifted_soliton = working_directory + "/espacio_parametros/mov_400cu/a=12.16_f=13.30_1"
    response_dir = working_directory + "/sigma_14mm"

    Z_shifted = np.loadtxt(shifted_soliton + '/Z_strobo.txt', delimiter=',')
    X_shifted = np.loadtxt(shifted_soliton + '/X_mm.txt', delimiter=',')
    T_shifted = np.loadtxt(shifted_soliton + '/T_strobo.txt', delimiter=',')
    X_mm = np.loadtxt(shifted_soliton + '/X_mm_image.txt', delimiter=',')
    img_shifted = plt.imread(shifted_soliton + '/image_modified.jpg')
    Z_shifted = filtro_superficie(Z_shifted, 8, "Y")
    Z_shifted = filtro_superficie(Z_shifted, 10, "X")

    response = np.loadtxt(response_dir + '/response.txt', delimiter=',')
    response_fit = np.loadtxt(response_dir + '/response_fit.txt', delimiter=',')
    response_X = np.loadtxt(response_dir + '/X_mm.txt', delimiter=',')

    img_shifted = cv2.cvtColor(img_shifted, cv2.COLOR_BGR2GRAY)

    Ny_img = len(img_shifted[:, 0])

    ti, tf = 0, 25
    xi, xf = -60, 60
    Ii, If = int(0.07 * Ny_img), int(0.47 * Ny_img)
    ji, jf = np.argmin(np.abs(X_mm - xi)), np.argmin(np.abs(X_mm - xf))

    fig, ax = plt.subplots(3, 1, figsize=(3, 3.5), gridspec_kw={'height_ratios': [1, 1, 1]}, dpi=300)
    ax[0].axis('off')
    ax[2].axis('off')

    # legend
    pcm = ax[1].pcolormesh(X_shifted, T_shifted, np.flip(Z_shifted, axis=1), cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm, aspect=10)
    cbar.set_label('$|A|$', rotation=0, size=15, labelpad=-14, y=1.3)
    cbar.set_ticks([0, 2, 4])

    ax[1].set_xlim([xi, xf])
    ax[1].tick_params(labelsize=15, direction='in', bottom=False, labelbottom=False)
    ax[1].set_ylim([0, 20])
    ax[1].set_ylabel('$t\ \\textrm{(s)}$', size=17)
    ax[1].grid(linestyle='--', alpha=0.2, color='k')
    cbar.ax.tick_params(labelsize=14)

    left, bottom, width, height = [0.231, 0.64, 0.576, 0.08]
    ax_up = fig.add_axes([left, bottom, width, height])
    ax_up.imshow(np.flip(img_shifted[Ii:If, ji:jf], axis=1), cmap="gray", aspect='auto')
    ax_up.tick_params(labelsize=15, left=False, right=False, top=False, bottom=False, labelleft=False, labelright=False, labeltop=False, labelbottom=False)

    left, bottom, width, height = [0.231, 0.28, 0.576, 0.08]
    ax_bot = fig.add_axes([left, bottom, width, height])
    ax_bot.plot(response_X, response, color='k', zorder=4)
    ax_bot.plot(response_X, response_fit, color='red', lw=0.8, zorder=5)
    ax_bot.hlines(0, xi, xf, colors="k", alpha=0.7, lw=0.5)
    ax_bot.set_xlim([xi, xf])
    ax_bot.set_xlabel(r'$x\ \textrm{(mm)}$', size=17)
    ax_bot.tick_params(labelsize=15, left=False, right=False, top=False, bottom=True, labelleft=False, labelright=False, labeltop=False, labelbottom=True)

    plt.tight_layout()
    plt.savefig("Fig01a.png", dpi=300)
    plt.close()