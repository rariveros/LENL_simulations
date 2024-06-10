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
    [alpha, beta, gamma_0, mu, nu, sigma_i, phi] = np.loadtxt(directory + '/parameters.txt', delimiter=',')
    # [alpha, beta, gamma_0, dist,  mu, nu, sigma]
    gamma_i = gamma_0 * (np.exp(-X ** 2 / (2 * sigma_i ** 2)))
    d = 4
    sigma = 10
    A = 0.2
    Z_r = Z_r / np.sqrt(beta)
    Z_i = Z_i / np.sqrt(beta)

    Z = Z_r + 1j * Z_i
    Z_conj = Z_r - 1j * Z_i
    Z_modulo = np.absolute(Z)

    phi_R = (A / 2) * (np.exp(- (X) ** 2 / (2 * sigma ** 2)) + np.exp(- (X) ** 2 / (2 * sigma ** 2))) * np.cos(np.sqrt(nu / alpha) * X)
    phi_I = (A / 2) * (np.exp(- (X + d/2) ** 2 / (2 * sigma ** 2)) - np.exp(- (X - d/2) ** 2 / (2 * sigma ** 2))) * np.sin(np.sqrt(nu / alpha) * X)
    plt.plot(X, Z_r, color="b")
    plt.plot(X, Z_i, color="r")
    plt.plot(X, phi_R, color="b", linestyle="--")
    plt.plot(X, phi_I, color="r", linestyle="--")
    #plt.plot(X, np.abs(Z), color="k")
    #plt.plot(X, gamma_i, color="k", linestyle="--")
    plt.grid(alpha=0.2)
    plt.show()