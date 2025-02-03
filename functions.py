from back_process import *

def equations_FD(eq, field_slices, t_i, x_grid, y_grid, parameters, operators):
    if eq == 'wave':
        U_1 = field_slices[0]
        U_2 = field_slices[1]
        DD = operators[0]
        c = parameters[0]

        ddU_x = Der(DD, U_1)
        ddU_y = np.transpose(Der(DD, np.transpose(U_1)))

        F = U_2
        G = c ** 2 * (ddU_x + ddU_y)

        fields = np.array([F, G])

    elif eq == 'wave_1D':
        U_1 = field_slices[0]
        U_2 = field_slices[1]
        DD = operators[0]
        c = parameters[0]

        ddU_x = Der(DD, U_1)

        F = U_2
        G = c ** 2 * (ddU_x)

        fields = np.array([F, G])

    elif eq == 'sine_gordon':
        U_1 = field_slices[0]
        U_2 = field_slices[1]
        DDx = operators[0]
        DDy = operators[1]
        c = parameters[0]
        gamma = parameters[1]
        alpha = parameters[2]

        #ddU_x = Der(DDx, U_1)
        ddU_x = np.transpose(Der(DDx, np.transpose(U_1)))

        F = U_2
        G = c ** 2 * (ddU_x) - gamma * U_2 + alpha * (U_1 - U_1 ** 3) / 3

        fields = np.array([F, G])

    elif eq == 'lug_lef':
        U_1 = field_slices[0]
        U_2 = field_slices[1]

        alpha = parameters[0]
        beta = parameters[1]
        S = parameters[2]
        S_1 = S[0]
        S_2 = S[1]
        mu = parameters[3]
        nu = parameters[4]

        DD = operators[0]

        ddU_1 = Der(DD, U_1)
        ddU_2 = Der(DD, U_2)

        F = S_1 - mu * U_1 + nu * U_2 + alpha * ddU_2 + beta * U_1 ** 2 * U_2 + U_2 ** 3
        G = S_2 - mu * U_2 - nu * U_1 - alpha * ddU_1 - beta * U_1 ** 3 - U_2 ** 2 * U_1

        fields = np.array([F, G])

    elif eq == 'KS':
        U_1 = field_slices[0]
        D1 = operators[0]
        D2 = operators[1]
        D4 = operators[2]

        dU_1 = Der(D1, U_1)
        ddU_1 = Der(D2, U_1)
        ddddU_1 = Der(D4, U_1)

        F = - ddddU_1 - ddU_1 - dU_1 ** 2

        fields = np.array([F])

    elif eq == 'PDNLS':

        U_1 = field_slices[0]
        U_2 = field_slices[1]

        alpha = parameters[0]
        beta = parameters[1]
        gamma = parameters[2]
        gamma_1 = gamma[0]
        gamma_2 = gamma[1]
        mu = parameters[3]
        nu = parameters[4]

        DD = operators[0]

        ddU_1 = Der(DD, U_1)
        ddU_2 = Der(DD, U_2)

        F = alpha * ddU_2 + (beta * (U_1 ** 2 + U_2 ** 2) + nu + gamma_2) * U_2 + (gamma_1 - mu) * U_1
        G = -alpha * ddU_1 - (beta * (U_1 ** 2 + U_2 ** 2) + nu - gamma_2) * U_1 - (gamma_1 + mu) * U_2

        fields = np.array([F, G])

    elif eq == 'choro_dinger':

        U_1 = field_slices[0]
        U_2 = field_slices[1]

        alpha = parameters[0]
        V = parameters[1]
        epsilon = 0.2
        V = V * (1 + epsilon * np.cos(2 * t_i))

        DD = operators[0]

        ddU_1 = Der(DD, U_1)
        ddU_2 = Der(DD, U_2)

        F = -alpha * ddU_2 + V * U_2
        G = alpha * ddU_1 - V * U_1

        fields = np.array([F, G])

    elif eq == 'PDNLS_complex':
        U = field_slices[0]

        alpha = parameters[0]
        beta = parameters[1]
        gamma = parameters[2]
        gamma = gamma[0]
        mu = parameters[3]
        nu = parameters[4]

        DD = operators[0]
        ddU = Der(DD, U)

        F = - (mu + 1j * nu) * U - 1j * alpha * ddU - 1j * beta * np.abs(U) * U + gamma * np.conjugate(U)

        fields = np.array([F])

    elif eq == 'PT_dimer':
        U1 = field_slices[0]
        U2 = field_slices[1]
        U1_conj = np.conjugate(U1)
        U2_conj = np.conjugate(U2)
        U1_mod = np.abs(U1)
        U2_mod = np.abs(U2)

        alpha = parameters[0]
        beta = parameters[1]
        mu = parameters[2]
        nu = parameters[3]
        sigma = parameters[4]
        pi = parameters[5]
        delta = parameters[6]

        S11 = sigma[0]
        S12 = sigma[1]
        S21 = sigma[2]
        S22 = sigma[3]

        P11 = pi[0]
        P12 = pi[1]
        P21 = pi[2]
        P22 = pi[3]

        D11 = delta[0]
        D21 = delta[1]
        D31 = delta[2]
        D41 = delta[3]
        D51 = delta[4]
        D61 = delta[5]

        D12 = delta[6]
        D22 = delta[7]
        D32 = delta[8]
        D42 = delta[9]
        D52 = delta[10]
        D62 = delta[11]

        F = - (mu + 1j * nu) * U1 - 1j * (S11 * U1 + S12 * U2) + P11 * U1_conj + P12 * U2_conj - 1j * beta * \
            (D11 * (U2_mod ** 2) * U2
             + D21 * (U2_mod ** 2) * U1
             + D31 * (U1_mod ** 2) * U2
             + D41 * (U1_mod ** 2) * U1
             + D51 * (U2 ** 2) * U1_conj
             + D61 * (U1 ** 2) * U2_conj
             )
        G = - (mu + 1j * nu) * U2 - 1j * (S21 * U1 + S22 * U2) + P21 * U1_conj + P22 * U2_conj - 1j * beta * \
            (D12 * (U2_mod ** 2) * U2
             + D22 * (U2_mod ** 2) * U1
             + D32 * (U1_mod ** 2) * U2
             + D42 * (U1_mod ** 2) * U1
             + D52 * (U2 ** 2) * U1_conj
             + D62 * (U1 ** 2) * U2_conj
             )
        fields = np.array([F, G])

    elif eq == 'soliton_variational':
        X = field_slices[0]
        Y = field_slices[1]

        alpha = parameters[0]
        beta = parameters[1]
        gamma_0 = parameters[2]
        mu = parameters[3]
        nu = parameters[4]
        sigma = parameters[5]
        dx = parameters[6]

        gammaR_ij = gamma_0 * (np.exp(-X ** 2 / (2 * sigma ** 2)))
        gammaI_ij = gamma_0 * (np.exp(-Y ** 2 / (2 * sigma ** 2)))
        phaseR_ij = 0.5 * np.arccos(mu / gammaR_ij)
        phaseI_ij = 0.5 * np.arccos(mu / gammaI_ij)
        deltaR_ij = 1 * (- nu + np.sqrt(gammaR_ij ** 2 - mu ** 2))
        deltaI_ij = 1 * (- nu + np.sqrt(gammaI_ij ** 2 - mu ** 2))
        AR_i = np.sqrt(2 * deltaR_ij) * np.cos(phaseR_ij)
        AI_i = - np.sqrt(2 * deltaI_ij) * np.sin(phaseI_ij)
        phi_01_i = (AR_i / np.cosh(np.sqrt(deltaR_ij / alpha) * (x_grid - X)))
        phi_02_i = (AI_i / np.cosh(np.sqrt(deltaI_ij / alpha) * (x_grid - Y)))
        Dphi_01_i = np.append(np.diff(phi_01_i) / dx, 0)
        Dphi_02_i = np.append(np.diff(phi_02_i) / dx, 0)
        f_ij = mu * phi_01_i * phi_02_i \
               - (- (alpha / 2) * Dphi_02_i ** 2 \
                  + (beta / 4) * phi_02_i ** 4 \
                  + (nu / 2) * phi_02_i ** 2 \
                  + gamma_0 * (np.exp(-x_grid ** 2 / (2 * sigma ** 2))) * phi_01_i * phi_02_i
                  + (beta / 2) * phi_01_i ** 2 * phi_02_i ** 2)
        F_ij = -integrate.simpson(f_ij, x_grid)
        g_ij = mu * phi_01_i * phi_02_i \
               + (- (alpha / 2) * Dphi_01_i ** 2 \
                  + (beta / 4) * phi_01_i ** 4 \
                  + (nu / 2) * phi_01_i ** 2 \
                  + gamma_0 * (np.exp(-x_grid ** 2 / (2 * sigma ** 2))) * phi_01_i * phi_02_i
                  + (beta / 2) * phi_01_i ** 2 * phi_02_i ** 2)
        G_ij = -integrate.simpson(g_ij, x_grid)
        h_ij = phi_01_i * phi_02_i
        H_ij = integrate.simpson(h_ij)

    elif eq == 'pdnlS_nospace':
        U1 = field_slices[0]

        alpha = parameters[0]
        beta = parameters[1]
        gamma = parameters[2]
        mu = parameters[3]
        nu = parameters[4]

        F = -(mu + 1j * nu) * U1 - 1j * np.abs(U1) ** 2 * U1 + gamma * np.conjugate(U1)
        fields = np.array([F])

    elif eq == 'PDNLS_interaction':

        U_1 = field_slices[0]
        U_2 = field_slices[1]

        V_1 = field_slices[2]
        V_2 = field_slices[3]

        alpha_1 = parameters[0]
        beta_1 = parameters[1]
        gamma1 = parameters[2]
        gamma1_1 = gamma1[0]
        gamma1_2 = gamma1[1]
        mu_1 = parameters[3]
        nu_1 = parameters[4]

        alpha_2 = parameters[5]
        beta_2 = parameters[6]
        gamma2 = parameters[7]
        gamma2_1 = gamma2[0]
        gamma2_2 = gamma2[1]
        mu_2 = parameters[8]
        nu_2 = parameters[9]

        DD = operators[0]

        ddU_1 = Der(DD, U_1)
        ddU_2 = Der(DD, U_2)
        ddV_1 = Der(DD, V_1)
        ddV_2 = Der(DD, V_2)

        interaction = 1

        F_1 = alpha_1 * ddU_2 + (beta_1 * (U_1 ** 2 + U_2 ** 2) + nu_1 + gamma1_2) * U_2 + (gamma1_1 - mu_1) * U_1 + interaction * (V_1 ** 2 + V_2 ** 2) * U_2
        G_1 = -alpha_1 * ddU_1 - (beta_1 * (U_1 ** 2 + U_2 ** 2) + nu_1 - gamma1_2) * U_1 - (gamma1_1 + mu_1) * U_2 - interaction * (V_1 ** 2 + V_2 ** 2) * U_1

        F_2 = alpha_2 * ddV_2 + (beta_2 * (V_1 ** 2 + V_2 ** 2) + nu_2 + gamma2_2) * V_2 + (gamma2_1 - mu_2) * V_1 + interaction * (U_1 ** 2 + U_2 ** 2) * V_2
        G_2 = -alpha_2 * ddV_1 - (beta_2 * (V_1 ** 2 + V_2 ** 2) + nu_2 - gamma2_2) * V_1 - (gamma2_1 + mu_2) * V_2 - interaction * (U_1 ** 2 + U_2 ** 2) * V_1

        fields = np.array([F_1, G_1, F_2, G_2])

    elif eq == "pdS":
        U_1 = field_slices[0]
        U_2 = field_slices[1]

        alpha = parameters[0]
        gamma = parameters[1]
        mu = parameters[2]
        nu = parameters[3]

        DD = operators[0]

        ddU_1 = Der(DD, U_1)
        ddU_2 = Der(DD, U_2)

        F = -alpha * ddU_2 - nu * U_2 + (gamma - mu) * U_1
        G = alpha * ddU_1 - nu * U_1 - (gamma + mu) * U_2

        fields = np.array([F, G])

    elif eq == 'kdV':

        U = field_slices[0]

        mu = parameters[0]

        D = operators[0]
        DDD = operators[1]

        dU = Der(D, U)
        dddU = Der(DDD, U)

        F = - 6 * U * dU - mu * dddU

        fields = np.array([F])
    elif eq == 'pdnlS_2D':

        U_1 = field_slices[0]
        U_2 = field_slices[1]

        alpha = parameters[0]
        beta = parameters[1]
        gamma = parameters[2]
        gamma_1 = gamma[0]
        gamma_2 = gamma[1]
        mu = parameters[3]
        nu = parameters[4]

        DDx = operators[0]
        DDy = operators[1]

        ddUx_1 = np.transpose(Der(DDx, np.transpose(U_1)))
        ddUy_1 = Der(DDy, U_1)

        ddUx_2 = np.transpose(Der(DDx, np.transpose(U_2)))
        ddUy_2 = Der(DDy, U_2)

        F = alpha * (ddUx_2 + ddUy_2) + (beta * (U_1 ** 2 + U_2 ** 2) + nu + gamma_2) * U_2 + (gamma_1 - mu) * U_1

        G = -alpha * (ddUx_1 + ddUy_1) - (beta * (U_1 ** 2 + U_2 ** 2) + nu - gamma_2) * U_1 - (gamma_1 + mu) * U_2

        fields = np.array([F, G])
    elif eq == 'PNDLS_gaussian_Xdamp':

        U_1 = field_slices[0]
        U_2 = field_slices[1]

        alpha = parameters[0]
        beta = parameters[1]
        gamma = parameters[2]
        gamma_1 = gamma[0]
        gamma_2 = gamma[1]
        mu_01 = parameters[3]
        mu_02 = parameters[4]
        nu = parameters[5]

        DD = operators[0]

        ddU_1 = Der(DD, U_1)
        ddU_2 = Der(DD, U_2)

        F = alpha * (ddU_2 - mu_02 * ddU_1) + (beta * (U_1 ** 2 + U_2 ** 2) + nu + gamma_2) * U_2 + (gamma_1 - mu_01) * U_1

        G = -alpha * (ddU_1 + mu_02 * ddU_2) - (beta * (U_1 ** 2 + U_2 ** 2) + nu - gamma_2) * U_1 - (gamma_1 + mu_01) * U_2

        fields = np.array([F, G])
    elif eq == 'pdnls_double_gaussian':
        U_1 = field_slices[0]
        U_2 = field_slices[1]
        V_1 = field_slices[2]
        V_2 = field_slices[3]

        alpha = parameters[0]
        beta = parameters[1]
        gamma = parameters[2]
        gamma_1 = gamma[0]
        gamma_2 = gamma[1]
        mu = parameters[3]
        nu = parameters[4]

        DD = operators[0]

        ddU_1 = Der(DD, U_1)
        ddU_2 = Der(DD, U_2)
        ddV_1 = Der(DD, V_1)
        ddV_2 = Der(DD, V_2)

        F_1 = alpha * ddU_2 + (beta * (U_1 ** 2 + U_2 ** 2) + nu + gamma_2) * U_2 + (gamma_1 - mu) * U_1 + 0.2 * (V_1 ** 2 + V_2 ** 2) * U_2

        G_1 = -alpha * ddU_1 - (beta * (U_1 ** 2 + U_2 ** 2) + nu - gamma_2) * U_1 - (gamma_1 + mu) * U_2 - 0.2 * (V_1 ** 2 + V_2 ** 2) * U_1

        F_2 = alpha * ddV_2 + (beta * (V_1 ** 2 + V_2 ** 2) + (nu + 0.1) + gamma_2) * V_2 + (gamma_1 - mu) * V_1 + 0.2 * (U_1 ** 2 + U_2 ** 2) * V_2

        G_2 = -alpha * ddV_1 - (beta * (V_1 ** 2 + V_2 ** 2) + (nu + 0.1) - gamma_2) * V_1 - (gamma_1 + mu) * V_2 - 0.2 * (U_1 ** 2 + U_2 ** 2) * V_1

        fields = np.array([F_1, G_1, F_2, G_2])
    elif eq == 'PNDLS_envelope_nonresonant':
        U_1 = field_slices[0]
        U_2 = field_slices[1]

        alpha = parameters[0]
        beta = parameters[1]
        gamma = parameters[2]
        mu = parameters[3]
        nu = parameters[4]
        sigma = parameters[5]

        delta = gamma - (mu + (1 / sigma) * np.sqrt(nu / alpha))
        mod_U = np.sqrt(U_1 ** 2 + U_2 ** 2)

        operator_1 = operators[0]
        operator_2 = operators[1]
        dU_1 = Der(operator_1, U_1)
        dU_2 = Der(operator_1, U_2)
        dUmod2 = Der(operator_1, mod_U ** 2)
        ddU_1 = Der(operator_2, U_1)
        ddU_2 = Der(operator_2, U_2)
        dU_1_mod = Der(operator_1, mod_U ** 2 * U_1)
        dU_2_mod = Der(operator_1, mod_U ** 2 * U_2)

        A = ((2 * nu) / (alpha * mu))
        B = ((delta - (mu / (2 * sigma ** 2)) * x_grid ** 2) - (9 * beta ** 2 * mod_U ** 4 / (2 * mu)))
        C = (2 * beta / mu) * np.sqrt(nu / alpha)

        F = A * ddU_1 + B * U_1 + C * (3 * dU_2_mod - 2 * U_2 * dUmod2) + (3 / (2 * np.sqrt(nu))) * \
            (np.cos(2 * np.sqrt(nu) * x_grid) * ((U_1 ** 2 - U_2 ** 2) * dU_2 + 2 * U_1 * U_2 * dU_1) +
             np.sin(2 * np.sqrt(nu) * x_grid) * (-1 * (U_1 ** 2 - U_2 ** 2) * dU_1 + 2 * U_1 * U_2 * dU_2))

        G = A * ddU_2 + B * U_2 - C * (3 * dU_1_mod - 2 * U_1 * dUmod2) - (3 * beta / (2 * np.sqrt(nu))) * \
            (np.cos(2 * np.sqrt(nu) * x_grid) * ((U_1 ** 2 - U_2 ** 2) * dU_1 - 2 * U_1 * U_2 * dU_2) +
             np.sin(2 * np.sqrt(nu) * x_grid) * ((U_1 ** 2 - U_2 ** 2) * dU_2 + 2 * U_1 * U_2 * dU_1))

        fields = np.array([F, G])

    elif eq == 'PNDLS_C_envelope':
        U_1 = field_slices[0]
        U_2 = field_slices[1]

        alpha = parameters[0]
        beta = parameters[1]
        gamma = parameters[2]
        mu = parameters[3]
        nu = parameters[4]
        sigma = parameters[5]

        delta = gamma - (mu + (1 / sigma) * np.sqrt(nu / alpha))
        mod_U = np.sqrt(U_1 ** 2 + U_2 ** 2)

        operator_1 = operators[0]
        operator_2 = operators[1]
        dUmod2 = Der(operator_1, mod_U ** 2)
        ddU_1 = Der(operator_2, U_1)
        ddU_2 = Der(operator_2, U_2)
        dU_1_mod = Der(operator_1, mod_U ** 2 * U_1)
        dU_2_mod = Der(operator_1, mod_U ** 2 * U_2)

        A = ((2 * nu) / (alpha * mu))
        B = ((delta - (mu / (2 * sigma ** 2)) * x_grid ** 2) - (9 * beta ** 2 * mod_U ** 4 / (2 * mu)))
        C = (2 * beta / mu) * np.sqrt(nu / alpha)

        F = A * ddU_1 + B * U_1 + C * (3 * dU_2_mod - 2 * U_2 * dUmod2)

        G = A * ddU_2 + B * U_2 - C * (3 * dU_1_mod - 2 * U_1 * dUmod2)

        fields = np.array([F, G])

    elif eq == 'LL_Zani':
        m1 = field_slices[0]
        m2 = field_slices[1]
        m3 = field_slices[2]

        alpha = parameters[0]
        C_ani = parameters[1]

        mod_23 = np.sqrt(m2 ** 2 + m3 ** 2)
        mod_13 = np.sqrt(m1 ** 2 + m3 ** 2)
        mod_12 = np.sqrt(m1 ** 2 + m2 ** 2)

        DD = operators[0]

        ddm1_x = Der(DD, m1)
        ddm1_y = np.transpose(Dxx(DD, np.transpose(m1)))
        ddm2_x = Der(DD, m2)
        ddm2_y = np.transpose(Dxx(DD, np.transpose(m2)))
        ddm3_x = Der(DD, m3)
        ddm3_y = np.transpose(Dxx(DD, np.transpose(m3)))

        h1 = ddm1_x + ddm1_y
        h2 = ddm2_x + ddm2_y
        h3 = ddm3_x + ddm3_y + C_ani * m3

        F = - (m2 * h3 - m3 * h2) + alpha * (mod_23 * h1 - m1 * m2 * h2 - m1 * m3 * h3)
        G = - (m3 * h1 - m1 * h3) + alpha * (mod_13 * h2 - m2 * m3 * h3 - m2 * m1 * h1)
        H = - (m1 * h2 - m2 * h1) + alpha * (mod_12 * h3 - m3 * m1 * h1 - m3 * m2 * h2)

        fields = np.array([F, G, H])

    elif eq == 'LL_x':
        m1 = field_slices[0]
        m2 = field_slices[1]
        m3 = field_slices[2]

        alpha = parameters[0]
        A = parameters[1]
        k_1 = parameters[2]
        k_2 = parameters[3]
        h_1 = parameters[4]
        h_2 = parameters[5]

        mod_23 = np.sqrt(m2 ** 2 + m3 ** 2)
        mod_13 = np.sqrt(m1 ** 2 + m3 ** 2)
        mod_12 = np.sqrt(m1 ** 2 + m2 ** 2)

        DD = operators[0]

        ddm1_x = Der(DD, m1)
        ddm2_x = Der(DD, m2)
        ddm3_x = Der(DD, m3)

        h1 = A * ddm1_x + k_1 * m1 + h_1
        h2 = A * ddm2_x - k_2 * m2 + h_2
        h3 = A * ddm3_x

        F = - (m2 * h3 - m3 * h2) + alpha * (mod_23 * h1 - m1 * m2 * h2 - m1 * m3 * h3)
        G = - (m3 * h1 - m1 * h3) + alpha * (mod_13 * h2 - m2 * m3 * h3 - m2 * m1 * h1)
        H = - (m1 * h2 - m2 * h1) + alpha * (mod_12 * h3 - m3 * m1 * h1 - m3 * m2 * h2)

    elif eq == 'LLGZ_01':
        m1 = field_slices[0]
        m2 = field_slices[1]
        m3 = field_slices[2]

        alpha = parameters[0]
        g = parameters[1]
        h_1 = parameters[2]
        h_d = parameters[3]
        A = parameters[4]

        mod_23 = np.sqrt(m2 ** 2 + m3 ** 2)
        mod_13 = np.sqrt(m1 ** 2 + m3 ** 2)
        mod_12 = np.sqrt(m1 ** 2 + m2 ** 2)

        DD = operators[0]

        ddm1_x = Der(DD, m1)
        ddm2_x = Der(DD, m2)
        ddm3_x = Der(DD, m3)

        h1 = A * ddm1_x + h_1
        h2 = A * ddm2_x
        h3 = A * ddm3_x - m3 * h_d

        F = - (m2 * h3 - m3 * h2) + alpha * (mod_23 * h1 - m1 * m2 * h2 - m1 * m3 * h3) - g * mod_23
        G = - (m3 * h1 - m1 * h3) + alpha * (mod_13 * h2 - m2 * m3 * h3 - m2 * m1 * h1) + g * m2 * m1
        H = - (m1 * h2 - m2 * h1) + alpha * (mod_12 * h3 - m3 * m1 * h1 - m3 * m2 * h2) + g * m3 * m1

        fields = np.array([F, G, H])

    elif eq == 'mathieu_single':
        U = field_slices[0]
        V = field_slices[1]

        w_0 = parameters[0]
        w_i = parameters[1]
        gamma_i = parameters[2]
        delta = parameters[3]

        F = V
        G = - (w_0 ** 2 + (gamma_i * w_i ** 2 / 2) * (np.sin(w_i * t_i) + np.sin((1 + delta) * w_i * t_i))) * (U - U ** 3 / 6)

        fields = np.array([F, G])

    elif eq == 'mathieu_double':
        U = field_slices[0]
        V = field_slices[1]

        w_0 = parameters[0]
        w_i = parameters[1]
        gamma_i = parameters[2]
        delta = parameters[3]

        interact = operators[0]
        OU = Der(interact, U)
        F = V
        G = - w_0 ** 2 * (U - U ** 3 / 6) + gamma_i * np.cos(w_i * t_i) + delta * OU

        fields = np.array([F, G])

    elif eq == "GPE":
        U_1 = field_slices[0]

        alpha = parameters[0]
        beta = parameters[1]
        V = parameters[2]
        mu = parameters[3]

        DD = operators[0]

        ddU_1 = Der(DD, U_1)

        F = -1j * (- alpha * ddU_1 + (beta * np.abs(U_1) ** 2 + V) * U_1)

        fields = np.array([F])

    elif eq == 'ladder_PDQHO':
        U1 = field_slices[0]
        U2 = field_slices[1]
        V1 = field_slices[2]
        V2 = field_slices[3]

        Delta = parameters[0]
        gamma = parameters[1]
        Omega = parameters[2]
        k = parameters[3]
        g = parameters[4]

        F1 = - (1j * Delta + gamma / 2) * U1 - 2 * 1j * Omega * U2 + 1j * k * V1 - 2 * 1j * g * np.abs(U1) ** 2 * U1
        F2 = + (1j * Delta - gamma / 2) * U2 + 2 * 1j * Omega * U1 - 1j * k * V2 + 2 * 1j * g * np.abs(U2) ** 2 * U2
        G1 = - (1j * Delta + gamma / 2) * V1 + 2 * 1j * Omega * V2 + 1j * k * U1 - 2 * 1j * g * np.abs(V1) ** 2 * V1
        G2 = + (1j * Delta - gamma / 2) * V2 - 2 * 1j * Omega * V1 - 1j * k * U2 + 2 * 1j * g * np.abs(V2) ** 2 * V2
        fields = np.array([F1, F2, G1, G2])

    elif eq == 'reduced_soliton':
        U = field_slices[0]
        V = field_slices[1]

        [alpha, beta, mu, nu, sigma, gamma] = parameters
        [D1, D2] = operators
        F1 = D1[0] + D1[1] * U + D1[2] * V + (1 / 2) * (D1[3] * U ** 2 + 2 * D1[4] * U * V + D1[5] * V ** 2) + (1 / 6) * (D1[6] * U ** 3 + 3 * D1[7] * U ** 2 * V + 3 * D1[8] * U * V ** 2 + D1[9] * V ** 3)
        F2 = D2[0] + D2[1] * U + D2[2] * V + (1 / 2) * (D2[3] * U ** 2 + 2 * D2[4] * U * V + D2[5] * V ** 2) + (1 / 6) * (D2[6] * U ** 3 + 3 * D2[7] * U ** 2 * V + 3 * D2[8] * U * V ** 2 + D2[9] * V ** 3)
        fields = np.array([F1, F2])
    return fields


def equations_FFT(eq, field_slices, parameters, x_grid, kappa):
    if eq == 'PNDLS_forced':
        U_1 = field_slices[0]
        U_2 = field_slices[1]

        alpha = parameters[0]
        beta = parameters[1]
        gamma = parameters[2]
        gamma_1 = gamma[0]
        gamma_2 = gamma[1]
        mu = parameters[3]
        nu = parameters[4]

        Uhat_1 = np.fft.fft(U_1)
        Uhat_2 = np.fft.fft(U_2)
        dd_Uhat_1 = -np.power(kappa, 2) * Uhat_1
        dd_Uhat_2 = -np.power(kappa, 2) * Uhat_2
        U_1 = np.fft.ifft(Uhat_1)
        U_2 = np.fft.ifft(Uhat_2)
        ddU_1 = np.fft.ifft(dd_Uhat_1)
        ddU_2 = np.fft.ifft(dd_Uhat_2)

        F = alpha * ddU_2 + (beta * (U_1 ** 2 + U_2 ** 2) + nu + gamma_2) * U_2 + (gamma_1 - mu) * U_1
        G = -alpha * ddU_1 - (beta * (U_1 ** 2 + U_2 ** 2) + nu - gamma_2) * U_1 - (gamma_1 + mu) * U_2
    elif eq == 'PDNLS':
        U_1 = field_slices[0]
        U_2 = field_slices[1]

        alpha = parameters[0]
        beta = parameters[1]
        gamma_0 = parameters[2]
        mu = parameters[3]
        nu = parameters[4]

        Uhat_1 = np.fft.fft(U_1)
        Uhat_2 = np.fft.fft(U_2)
        dd_Uhat_1 = -np.power(kappa, 2) * Uhat_1
        dd_Uhat_2 = -np.power(kappa, 2) * Uhat_2
        U_1 = np.fft.ifft(Uhat_1)
        U_2 = np.fft.ifft(Uhat_2)
        ddU_1 = np.fft.ifft(dd_Uhat_1)
        ddU_2 = np.fft.ifft(dd_Uhat_2)

        F = alpha * ddU_2 + (beta * (U_1 ** 2 + U_2 ** 2) + nu) * U_2 + (gamma_0 - mu) * U_1
        G = - alpha * ddU_1 - (beta * (U_1 ** 2 + U_2 ** 2) + nu) * U_1 - (gamma_0 + mu) * U_2

    elif eq == 'PDNLS_drift':
        U_1 = field_slices[0]
        U_2 = field_slices[1]

        alpha = parameters[0]
        beta = parameters[1]
        gamma_0 = parameters[2]
        mu = parameters[3]
        nu = parameters[4]

        Uhat_1 = np.fft.fft(U_1)
        Uhat_2 = np.fft.fft(U_2)
        dd_Uhat_1 = -np.power(kappa, 2) * Uhat_1
        dd_Uhat_2 = -np.power(kappa, 2) * Uhat_2
        U_1 = np.fft.ifft(Uhat_1)
        U_2 = np.fft.ifft(Uhat_2)
        ddU_1 = np.fft.ifft(dd_Uhat_1)
        ddU_2 = np.fft.ifft(dd_Uhat_2)

        F = alpha * ddU_2 + (beta * (U_1 ** 2 + U_2 ** 2) + nu) * U_2 + (gamma_0 * (1 + 0.01 * (x_grid - x_grid[0])) - mu) * U_1
        G = - alpha * ddU_1 - (beta * (U_1 ** 2 + U_2 ** 2) + nu) * U_1 - (gamma_0 * (1 + 0.01 * (x_grid - x_grid[0])) + mu) * U_2
    elif eq == 'drift_example':
        U_1 = field_slices[0]
        U_2 = field_slices[1]

        alpha = parameters[0]
        beta = parameters[1]
        mu = parameters[2]
        v_0 = parameters[3]
        nu = parameters[4]

        mod_U = np.sqrt(U_1 ** 2 + U_2 ** 2)
        Uhat_1 = np.fft.fft(U_1)
        Uhat_2 = np.fft.fft(U_2)
        dd_Uhat_1 = -np.power(kappa, 2) * Uhat_1
        dd_Uhat_2 = -np.power(kappa, 2) * Uhat_2
        U_1 = np.fft.ifft(Uhat_1)
        U_2 = np.fft.ifft(Uhat_2)
        ddU_1 = np.fft.ifft(dd_Uhat_1)
        ddU_2 = np.fft.ifft(dd_Uhat_2)
        sigma = 10
        F = (mu - beta * mod_U) * U_1 + alpha * ddU_1 + (nu + v_0 * np.exp(- x_grid ** 2 / 2 * sigma) ) * U_2
        G = (mu - beta * mod_U) * U_2 + alpha * ddU_2 - (nu + v_0 * np.exp(- x_grid ** 2 / 2 * sigma) ) * U_1
    return np.array([F.real, G.real])


def equations_ode(eq, vect, t, parameters):
    if eq == 'LL_ode':
        mx = vect[0]
        my = vect[1]
        mz = vect[2]
        ms = np.sqrt(mx ** 2 + my ** 2 + mz ** 2)

        gamma_L = parameters[1]
        alpha = parameters[2]
        omega_x = parameters[3]
        omega_y = parameters[4]
        delta_x = parameters[5]
        delta_y = parameters[6]

        H = parameters[0]
        hx = H[0]
        hx = hx * delta_x * np.sin(omega_x * t)
        hy = H[1]
        hy = hy * delta_y * np.sin(omega_y * t)
        hz = H[2]

        a = gamma_L
        b = alpha * gamma_L / ms

        f = -a * (my * hz - mz * hy) + b * ((my ** 2 + mz ** 2) * hx - mx * my * hy - mx * mz * hz)
        g = -a * (mz * hx - mx * hz) + b * ((mx ** 2 + mz ** 2) * hy - my * mz * hz - my * mx * hx)
        h = -a * (mx * hy - my * hx) + b * ((mx ** 2 + my ** 2) * hz - mz * mx * hx - mz * my * hy)
        vects = np.array([f, g, h])
    elif eq == "transcritical":
        U_1 = vect[0]

        alpha = parameters[0]
        beta = parameters[1]
        gamma = parameters[2]

        noise = 2 * (np.random.rand(1) - 0.5)

        F = alpha * U_1 - U_1 ** 2 + np.sqrt(beta) * noise
        vects = np.array([F])
    return vects