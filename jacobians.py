from back_process import *

def jacobians_FD(eq, fields, t_grid, x_grid, y_grid, parameters, operators):
    if eq == 'pdnlS' or eq == 'PDNLS':
        U_1 = fields[0]
        U_2 = fields[1]
        [alpha, beta, gamma, mu, nu] = parameters
        #print(alpha)
        #print(beta)
        #print(nu)
        #print(mu)
        gamma_1 = gamma[0]

        DD = operators[0]

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

        J_1 = np.concatenate((J_11, J_12), axis=1)
        J_2 = np.concatenate((J_21, J_22), axis=1)
        J = np.concatenate((J_1, J_2), axis=0)

    elif eq == 'lug_lef':
        U_1 = fields[0]
        U_2 = fields[1]
        delta = parameters[0]
        DD = operators[0]

        j_11 = [-1 + 2 * U_1 * U_2]
        j_12 = [delta + U_1 ** 2 + 3 * U_2 ** 2]
        j_21 = [-delta - U_2 ** 2 - 3 * U_1 ** 2]
        j_22 = [-1 - 2 * U_2 * U_1]

        J_11 = diags(j_11, [0])
        J_11 = J_11.toarray()

        J_12 = diags(j_12, [0]) + DD
        J_12 = J_12.toarray()

        J_21 = diags(j_21, [0]) - DD
        J_21 = J_21.toarray()

        J_22 = diags(j_22, [0])
        J_22 = J_22.toarray()

        J_1 = np.concatenate((J_11, J_12), axis=1)
        J_2 = np.concatenate((J_21, J_22), axis=1)
        J = np.concatenate((J_1, J_2), axis=0)
    return J
