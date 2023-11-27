import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    disco = 'E:/'
    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FD'
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    directories = os.listdir(working_directory)
    arg_max = []
    mod_max = []
    x_arg_max = []
    x_mod_max = []
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    N_dir = len(directories)
    sigmas = np.arange(4, 25, 1)
    for i in range(N_dir):
        directory_i = directories[i]
        directory = working_directory + "/" + directory_i
        files = os.listdir(directory)
        Z_r = np.loadtxt(directory + '/gamma=0.150/field_real.txt', delimiter=',')
        Z_i = np.loadtxt(directory + '/gamma=0.150/field_img.txt', delimiter=',')
        X = np.loadtxt(directory + '/gamma=0.150/X.txt', delimiter=',')
        T = np.loadtxt(directory + '/gamma=0.150/T.txt', delimiter=',')
        params = np.loadtxt(directory + '/gamma=0.150/parameters.txt', delimiter=',')

        Z_complex = Z_r + 1j * Z_i
        Z_complex = Z_complex[:, 40:-40]
        X = X[40:-40]
        arg = np.angle(Z_complex)
        arg = (2 * np.pi + arg) * (arg < 0) + arg * (arg > 0)
        Z_modulo = np.absolute(Z_complex)
        arg = np.unwrap(arg[-1, :], period=2 * np.pi)

        l_cut = 80
        r_cut = -80
        arg_left = arg[:l_cut]
        arg_right = arg[r_cut:]
        x_left = X[:l_cut]
        x_right = X[r_cut:]
        slope_left = np.mean(np.gradient(arg_left, x_left))
        slope_std_left = np.std(np.gradient(arg_left, x_left))
        slope_right = np.mean(np.gradient(arg_right, x_right))
        slope_std_right = np.std(np.gradient(arg_right, x_right))

        n_color = (i / N_dir)
        ax1.plot(X, arg + 0.1 * i, c=(n_color, 0, 1-n_color), linewidth=1)
        ax2.plot(X, (Z_modulo[-1, :] / np.amax(Z_modulo[-1, :])) + 0.1 * i, label="$\sigma_i = " + str(np.round(np.abs(sigmas[i]), 3)) + "\ \\textrm{mm}$", c=(n_color, 0, 1-n_color), linewidth=1)
        ax3.errorbar(np.round(np.abs(sigmas[i]), 3), slope_left, slope_std_left, c=(n_color, 0, 1-n_color), fmt='o')
        ax4.errorbar(np.round(np.abs(sigmas[i]), 3), slope_right, slope_std_right, c=(n_color, 0, 1-n_color), fmt='o')

        arg_max.append(np.amax(arg))
        x_arg_max.append(X[np.argmax(arg)])
        mod_max.append(np.amax(Z_modulo[-1, :]))
        x_mod_max.append(X[np.argmax(Z_modulo[-1, :])])
    ti, tf = 720, 1000
    xi, xf = -100, 100

    lim_left = np.ones(len(x_left))
    lim_right = np.ones(len(x_right))
    ax1.set_xlim(xi, xf)
    ax1.set_ylim(0, 10)
    ax1.fill_between(x_left, 10 * lim_left, alpha=0.3, color="r", edgecolor="none")
    ax1.fill_between(x_right, 10 * lim_right, alpha=0.3, color="r", edgecolor="none")
    ax1.vlines([x_arg_max[0], x_arg_max[-1]], 0, 10, linestyles="--", colors="k", linewidths=1)
    ax1.text(0.1, 0.2, '$\\theta_{L}(x)$', horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes)
    ax1.text(0.9, 0.2, '$\\theta_{R}(x)$', horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes)
    ax1.text(0.35, 0.35, '$\delta x_{\\theta}= ' + str(np.abs(x_arg_max[0] - x_arg_max[-1])) + '$', horizontalalignment='center',
             verticalalignment='center', transform=ax1.transAxes)
    ax1.grid(alpha=0.4)
    ax1.set_ylabel("$\\theta(x)$", fontsize=12)
    ax1.set_xlabel("$x\ \\textrm{(mm)}$", fontsize=12)

    ax2.vlines([x_mod_max[0], x_mod_max[-1]], 0, 7, linestyles="--", colors="k", linewidths=1)
    ax2.legend(loc="upper right", fontsize=5)
    ax2.text(0.17, 0.85, '$\delta x_{|A|}= ' + str(np.abs(x_mod_max[0] - x_mod_max[-1])) + '$', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes)
    ax2.set_ylim(0, 4)
    ax2.set_xlim(xi/2, xf/2)
    ax2.grid(alpha=0.4)
    ax2.set_ylabel("$|A(x)|$", fontsize=12)
    ax2.set_xlabel("$x\ \\textrm{(mm)}$", fontsize=12)

    ax3.grid(alpha=0.4)
    ax3.set_ylabel("$\partial_x \\theta_L(x)$", fontsize=12)
    ax3.set_xlabel("$|\\nu|$", fontsize=12)
    ax3.set_ylim(0.028, 0.053)

    ax4.grid(alpha=0.4)
    ax4.set_ylabel("$\partial_x \\theta_R(x)$", fontsize=12)
    ax4.set_xlabel("$|\\nu|$", fontsize=12)
    ax4.set_ylim(-0.053, -0.028)
    plt.tight_layout()
    plt.savefig("preliminary_sigma.png", dpi=300)