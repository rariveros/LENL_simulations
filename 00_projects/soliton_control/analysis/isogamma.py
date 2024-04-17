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
    f_i_01 = np.arange(13.1, 14.3, 0.05)
    a_01 = 0.19 / (f_i_01 ** 2 * 8.401093 * 10 ** (-5))

    f_i_02 = np.arange(13.1, 14.3, 0.05)
    a_02 = 0.18 / (f_i_02 ** 2 * 8.401093 * 10 ** (-5))
    #alpha, beta, nu, gamma, f_0 = fluid_pdnls_parameters(f_i_02, a_02, 20)

    f_i_03 = np.arange(13.1, 14.3, 0.05)
    a_03 = 0.17 / (f_i_03 ** 2 * 8.401093 * 10 ** (-5))
    alpha, beta, nu, gamma, f_0 = fluid_pdnls_parameters(f_i_03, a_03, 20)
    print(np.flip(a_03))
    print(np.flip(f_i_03))