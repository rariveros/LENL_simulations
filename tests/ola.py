from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    data = np.loadtxt('C:/Users/Rafael Riveros/Desktop/ola.txt',dtype=str)
    data_numerical = []
    names = []
    for i in range(len(data[:, 0])):
        data_numerical.append(data[i, 1:])
        names.append(data[i, 0])
    data_numerical = np.array(data_numerical)
