from back_process import *

def fluid_pdnls_parameters(f_i, a_ang, d):
    g = 9790
    l_y = 14.2
    w = 2 * np.pi * (f_i / 2)
    k_y = np.pi / l_y
    k = k_y
    tau = np.tanh(k * d)
    w_1 = np.sqrt(g * k * tau)
    f_0 = w_1 / np.pi
    gamma = f_i ** 2 * a_ang * 8.401093 * 10 ** (-5)
    alpha = (1 / (4 * k ** 2)) * (1 + k * d * ((1 - tau ** 2) / tau))
    beta = (k ** 2 / 64) * (6 * tau ** 2 - 5 + 16 * tau ** (-2) - 9 * tau ** (-4))  # término no lineal
    nu = 0.5 * ((w / w_1) ** 2 - 1)
    return alpha, beta, nu, gamma, f_0

if __name__ == "__main__":
    a_0 = 12
    delta_x = np.array([10, 9, 8.5, 8, 7, 6.5, 6, 5, 3, 0])
    d = 20
    f_0 = 13.5
    alpha, beta, nu, gamma, f_0 = fluid_pdnls_parameters(f_0, a_0, d)
    f_i = np.arange(13.5, 14, 0.05)
    a = gamma / (f_i ** 2 * 8.401093 * 10 ** (-5))
    alpha, beta, nu, gamma, f_0 = fluid_pdnls_parameters(f_i, a, d)
    plt.scatter(nu, delta_x, c="k")
    plt.xlabel("$\\nu$", fontsize=20)
    plt.ylabel("$\Delta x$", fontsize=20)
    plt.xlim([-0.1, 0])
    #plt.ylim([0, 0.2])
    plt.grid(alpha=0.2)
    plt.show()

    plt.show()