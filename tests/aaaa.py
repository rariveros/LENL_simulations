import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    xd = np.arange(973600)
    xdd = np.arange(311000) + 311000
    a = []
    for i in range(len(xd)):
        if xd[i] > xdd[0] and xd[i] < xdd[-1]:
            a.append(xd[i])
    print(a)