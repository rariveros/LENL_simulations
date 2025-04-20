import matplotlib.pyplot as plt
import numpy as np

from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    gamma = np.arange(0.2, 0.35, 0.005)
    mu = 0.1
    nu = 0
    k = 0.1
    R0 = P0 = np.sqrt(-nu + np.sqrt(gamma ** 2 - mu ** 2))
    R1 = (R0 / gamma) * ((np.sqrt(gamma ** 2 - mu ** 2) + mu) / (- 2 * nu + 4 * np.sqrt(gamma ** 2 - mu ** 2)))
    P1 = -(P0 / gamma) * ((np.sqrt(gamma ** 2 - mu ** 2) + mu) / (- 2 * nu + 4 * np.sqrt(gamma ** 2 - mu ** 2)))
    R = R0 + 0.01 * k * R1  # 0.336 #
    P = P0 + 0.01 * k * P1  # 0.23 #

    theta = np.arctan(P / R)
    phi = np.arctan(R / P)
    thetab = np.arctan(-P / R)
    phib = np.arctan(-R / P)

    plt.scatter(gamma, theta, color='b')
    plt.scatter(gamma, phi, color='r')
    plt.scatter(gamma, thetab, color='b')
    plt.scatter(gamma, phib, color='r')
    plt.show()
    plt.close()

    u_abs2 = ((gamma + mu) / gamma) * (-nu + np.sqrt(gamma ** 2 - mu ** 2))
    v_abs2 = u_abs2 * (gamma - mu) / (gamma + mu)
    plt.plot(gamma, np.sqrt(u_abs2), color='b')
    plt.plot(gamma, np.sqrt(v_abs2), color='r')
    plt.show()