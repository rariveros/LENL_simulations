from back_process import *

def hamiltonians_FD(eq, fields, t_grid, x_grid, y_grid, parameters, operators):
    if eq == 'PDNLS':
        Z_r = fields[0]
        Z_i = fields[1]
        D = operators[1]
        alpha = parameters[0]
        beta = parameters[1]
        gamma = parameters[2]
        gamma_1 = gamma[0]
        gamma_2 = gamma[1]
        mu = parameters[3]
        nu = parameters[4]

        hamiltonian_x = Der(D, Z_r) ** 2 + Der(D, Z_i) ** 2 - 0.5 * (Z_r ** 2 + Z_i ** 2) ** 2 - nu * (
                    Z_r ** 2 + Z_i ** 2)
        hamiltonian = integrate.simpson(hamiltonian_x, x_grid)
    return hamiltonian
