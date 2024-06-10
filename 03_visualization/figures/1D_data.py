from back_process import *


if __name__ == '__main__':
    disco = 'C:/'
    initial_dir_data = str(disco) + ':/mnustes_science/simulation_data'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    Z_r = np.loadtxt(directory + '/field_real_0.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img_0.txt', delimiter=',')
    X = np.loadtxt(directory + '/X.txt', delimiter=',')
    [alpha, beta, gamma_0, mu, nu, sigma] = np.loadtxt(directory + '/parameters.txt', delimiter=',')
    # [alpha, beta, gamma_0, dist,  mu, nu, sigma]
    gamma_i = gamma_0 * (np.exp(-X ** 2 / (2 * sigma ** 2)))

    Z_r = Z_r / np.sqrt(beta)
    Z_i = Z_i / np.sqrt(beta)

    Z = Z_r + 1j * Z_i
    Z_conj = Z_r - 1j * Z_i
    Z_modulo = np.absolute(Z)
    plt.plot(X, Z_r[-1, :],color="b")
    plt.plot(X, Z_i[-1, :],color="r")
    plt.plot(X, np.abs(Z[-1, :]), color="k")
    plt.plot(X, gamma_i, color="k", linestyle="--")
    plt.hlines(0.075, -120, 120, colors="k", linestyles="--")
    plt.hlines(0., -120, 120, colors="k")
    plt.grid(alpha=0.2)
    plt.show()