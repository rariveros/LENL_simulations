from back_process import *
from jacobians import *

if __name__ == '__main__':
    disco = 'aD:/'
    eq = 'pdnlS'

    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FD'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    x_grid = np.loadtxt(directory + '/X.txt', delimiter=',')
    t_grid = np.loadtxt(directory + '/T.txt', delimiter=',')
    [alpha, beta, gamma_0, dist, sigma,  mu, nu] = np.loadtxt(directory + '/parameters.txt', delimiter=',')
    print("[alpha, beta, gamma_0, mu, nu, sigma]")
    print([alpha, beta, gamma_0, mu, nu, sigma])
    Nx = len(x_grid)
    dx = x_grid[1] - x_grid[0]
    phi = np.pi
    gamma_complex = gamma_0 * (np.exp(- (x_grid - dist / 2) ** 2 / (2 * sigma ** 2)) + np.exp(1j * phi) * np.exp(- (x_grid + dist / 2) ** 2 / (2 * sigma ** 2)))
    gamma_real = np.real(gamma_complex)
    gamma_img = np.imag(gamma_complex)
    gamma = np.array([gamma_real, gamma_img])
    parameters = [alpha, beta, gamma, mu, nu]
    Z = Z_r + 1j * Z_i
    Z_conj = np.conj(Z)
    Ks = np.arange(-2, 2, 0.002)
    Zk = []
    for k in Ks:
        Zk_i = (1 / (2 * np.pi)) * integrate.simpson(Z[-1, :] * np.exp(- 1j * k * x_grid), x_grid)
        Zk.append(Zk_i)
    plt.plot(Ks, np.abs(Zk), color="k")
    plt.show()

    plt.plot(x_grid, np.real(gamma_complex), color="k")
    plt.plot(x_grid, Z_r[-1, :], color="b")
    plt.plot(x_grid, Z_i[-1, :], color="r")
    plt.show()