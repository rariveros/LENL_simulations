from functions import *
from back_process import *
from time_integrators import *

def gaussian(x, A, x_0, sigma):
    return A * np.exp(- (x - x_0) ** 2 / (2 * sigma ** 2))

if __name__ == '__main__':
    x_grid = np.arange(-75, 75, 0.1)
    gamma0 = 0.2
    mu = 0.1
    d = 20
    sigma = 15
    gamma_L = gaussian(x_grid, gamma0, -d / 2, sigma)
    gamma_R = -gaussian(x_grid, gamma0, d / 2, sigma)
    gamma = gamma_L + gamma_R
    plt.plot(x_grid, gamma - mu)
    plt.plot(x_grid, -(gamma + mu))
    plt.show()