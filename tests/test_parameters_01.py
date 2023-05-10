from functions import *
from back_process import *
from time_integrators import *
from directories_lyap import *

if __name__ == '__main__':
    f_i = 15.195
    a = 14.5
    d = 20
    alpha, beta, nu, gamma = fluid_pdnls_parameters(f_i, a, d)
    print(alpha, beta, nu, gamma)