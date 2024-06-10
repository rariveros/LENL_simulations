from back_process import *


if __name__ == '__main__':
    directory = 'C:/mnustes_science/simulation_data/FD/variational_pdnlS/localized_patterns/alpha=6.524/beta=1.000/mu=0.100/nu=0.018/sigma=6.000/gamma=0.185'

    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    X = np.loadtxt(directory + '/X.txt', delimiter=',')
    params = np.loadtxt(directory + '/parameters.txt', delimiter=',')
    # [alpha, beta, gamma_0, dist,  mu, nu, sigma]

    def zeta_01(x, A, sigma, k):
        return A * (np.exp((- x ** 2) / (2 * sigma ** 2))) * np.cos(k * x)

    def zeta_02(x, A, sigma, k):
        return -A * np.sqrt(1 / (2 * sigma ** 2)) * x * (np.exp((- x ** 2) / (2 * sigma ** 2))) * np.sin(k * x)

    param_01, cov_01 = curve_fit(zeta_01, X, Z_r, bounds=[(0.1, 0, 0.01), (0.5, 25, 0.2)])
    param_02, cov_02 = curve_fit(zeta_02, X, Z_i, bounds=[(0.1, 0, 0.01), (0.5, 25, 0.2)])


    A_01 = param_01[0]
    sigma_01 = param_01[1]
    k_01 = param_01[2]

    A_02 = param_02[0]
    sigma_02 = param_02[1]
    k_02 = param_02[2]
    print(A_01, sigma_01, k_01)
    print(A_02, sigma_02, k_02)
    fit_01 = zeta_01(X, A_01, sigma_01, k_01)
    #fit_02 = zeta_02(X, A_02, sigma_02, k_02)
    #fit_01 = zeta_01(X, A_02, x0_02, sigma_02, k_02)
    fit_02 = zeta_02(X, A_01, sigma_01, k_01)
    plt.plot(X, Z_r, color="b")
    plt.plot(X, Z_i, color='r')
    plt.plot(X, fit_01, '--', color="b")
    plt.plot(X, fit_02, '--', color='r')
    plt.show()