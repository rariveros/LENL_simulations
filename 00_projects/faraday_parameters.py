from functions import *
from back_process import *
from time_integrators import *
from directories_lyap import *

if __name__ == '__main__':
    f_i = 14.80
    a = 11.4
    d = 20
    alpha, beta, nu, gamma = fluid_pdnls_parameters(f_i, a, d)
    print(alpha, beta, nu, gamma)