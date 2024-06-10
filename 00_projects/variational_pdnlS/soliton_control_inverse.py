from back_process import *


if __name__ == '__main__':
    disco = 'C:/'
    initial_dir_data = str(disco) + ':/mnustes_science/simulation_data'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    x_grid = np.loadtxt(directory + '/X.txt', delimiter=',')
    T = np.loadtxt(directory + '/T.txt', delimiter=',')
    [alpha, beta, gamma_0, mu, nu, sigma_i] = np.loadtxt(directory + '/parameters.txt', delimiter=',')
    # [alpha, beta, gamma_0, mu, nu, sigma]

    Z_r = Z_r[-1, :]
    Z_i = Z_i[-1, :]

    dx = x_grid[1] - x_grid[0]

    X = np.arange(-20, 20, 0.1)

    def soliton(x, A, b, x0):
        return A * (1 / np.cosh(b * (x - x0)))


    p1, c1 = curve_fit(soliton, x_grid, Z_r)
    p2, c2 = curve_fit(soliton, x_grid, Z_i)
    A01 = p1[0]
    b01 = p1[1]
    x01 = p1[2]
    A02 = p2[0]
    b02 = p2[1]
    x02 = p2[2]
    F = []
    for i in range(len(X)):
        phi_01_i = soliton(x_grid, A01, b01, x01)
        phi_02_i = soliton(x_grid, A02, b02, x02)
        #plt.plot(x_grid, Z_r, color="b")
        #plt.plot(x_grid, Z_i, color="r")
        #plt.plot(x_grid, phi_01_i, color="b", linestyle="--")
        #plt.plot(x_grid, phi_02_i, color="r", linestyle="--")
        #plt.show()
        Dphi_01_i = np.append(np.diff(phi_01_i) / (dx), 0)
        Dphi_02_i = np.append(np.diff(phi_02_i) / (dx), 0)
        f_i = (-(alpha / 2) * (Dphi_01_i ** 2 + Dphi_02_i ** 2) + (beta / 4) * (phi_01_i ** 2 + phi_02_i ** 2) ** 2 + (nu / 2) * (phi_01_i ** 2 + phi_02_i ** 2) + gamma_0 * (np.exp(-(x_grid - X[i]) ** 2/(2 * sigma_i ** 2))) * phi_01_i * phi_02_i) / (2 * mu)
        F_i = integrate.simpson(f_i, x_grid)
        F.append(F_i)
    F = np.array(F)
    dF = np.diff(F) / (X[1] - X[0])
    ddF = np.diff(dF) / (X[1] - X[0])
    argmin = np.argmax(np.real(F))
    print(X[argmin])
    plt.plot(X, np.real(F))
    #plt.plot(X, np.imag(ddF))
    plt.show()
    #X, Y = np.meshgrid(b, A)
    #fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    #surf = ax.plot_surface(X, Y, F, cmap=parula_map, linewidth=0, antialiased=False)
    #ax.set_xlabel("$b$")
    #ax.set_ylabel("$A$")
    #plt.show()