import matplotlib.pyplot as plt
from time_integrators import *

def gaussian(x, A, x_0, sigma):
    return A * np.exp(- (x - x_0) ** 2 / (2 * sigma ** 2))

if __name__ == '__main__':
    eq = 'pdnlS'
    [alpha, beta, gamma_0, mu, nu, sigma] = [6.524, 1, 0.185, 0.1, 0.018, 6]
    dx = 1
    nu = 0
    [xmin, xmax, dx] = [-50, 50, dx]
    x_grid = np.arange(xmin, xmax, dx)
    dist = 25
    Nx = len(x_grid)
    dx = x_grid[1] - x_grid[0]
    gamma_l = 1
    gamma_r = 1
    gamma_img = 0
    gamma_real = -(gaussian(x_grid, gamma_0, -dist / 2, sigma) - gaussian(x_grid, gamma_0, dist / 2, sigma))
    gamma = [gamma_real, gamma_img]
    parameters = [alpha, beta, gamma, mu, nu]
    DD = sparse_DD_neumann(Nx, dx)

    [alpha, beta, gamma, mu, nu] = parameters
    # print(alpha)
    # print(beta)
    # print(nu)
    # print(mu)
    alpha = 1
    gamma_1 = gamma_real
    U_1 = 0 * x_grid
    U_2 = 0 * x_grid

    j_11 = [gamma_1 - mu + 2 * beta * U_1 * U_2]
    j_12 = [nu + beta * U_1 ** 2 + 3 * beta * U_2 ** 2]
    j_21 = [- nu - beta * U_2 ** 2 - 3 * beta * U_1 ** 2]
    j_22 = [- (gamma_1 + mu) - 2 * beta * U_2 * U_1]

    J_11 = diags(j_11, [0])
    J_11 = J_11.toarray()

    J_12 = diags(j_12, [0]) + alpha * DD
    J_12 = J_12.toarray()

    J_21 = diags(j_21, [0]) - alpha * DD
    J_21 = J_21.toarray()

    J_22 = diags(j_22, [0])
    J_22 = J_22.toarray()

    plt.imshow(J_12 @ J_21)
    plt.show()
    plt.close()
