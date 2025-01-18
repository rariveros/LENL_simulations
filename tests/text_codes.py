import matplotlib.pyplot as plt
from functions import *
from time_integrators import *
from back_process import *
import numpy as np

if __name__ == '__main__':
    a = np.arange(0, 10, 0.01)
    b = np.exp(- a ** 2)
    plt.plot(a, b)
    plt.show()