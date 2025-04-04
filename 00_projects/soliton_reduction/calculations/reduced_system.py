from functions import *
from back_process import *
from time_integrators import *


if __name__ == '__main__':
    dx = 0.1
    x_grid = np.arange(-40, 40, dx)
    Xs = np.arange(-0.1, 0.1, 0.02)
    Ys = np.arange(-0.2, 0.2, 0.02)
    dX = np.abs(Xs[1] - Xs[0])
    dY = np.abs(Ys[1] - Ys[0])
    NX = len(Xs)
    NY = len(Ys)
    I = np.argmin(np.abs(Xs))
    J = np.argmin(np.abs(Ys))
    print(Xs[I], Ys[J])
    nus = np.arange(-0.2, 0.0, 0.005)
    for nu in nus:
        print("nu =" + f"{nu:.{3}f}")
        [alpha, beta, mu, nu, sigma, gamma] = [6.524, 1, 0.075, nu, 15, 0.18]
        gamma_x = gamma * np.exp(- x_grid ** 2 / (2 * sigma ** 2))
        F_01 = []
        F_02 = []
        I1 = []
        for X in Xs:
            F_01_i = []
            F_02_i = []
            I1_i = []
            for Y in Ys:
                [phi_01, phi_02, Dphi_01, Dphi_02] = phis(alpha, beta, nu, mu, gamma, sigma, X, Y, x_grid, dx)
                #plt.plot(x_grid, phi_01, color="b")
                #plt.plot(x_grid, Dphi_01, color="b", linestyle="dashed")
                #plt.plot(x_grid, phi_02, color="r")
                #plt.plot(x_grid, Dphi_02, color="r", linestyle="dashed")
                #plt.show()
                i1 = 0.5 * phi_01 * phi_02
                i2a = 0.5 * (-alpha * Dphi_01 ** 2 + nu * phi_01 ** 2 + (beta / 2) * phi_01 ** 4)
                i2b = 0.5 * (-alpha * Dphi_02 ** 2 + nu * phi_02 ** 2 + (beta / 2) * phi_02 ** 4)
                i3 = (beta / 2) * phi_01 ** 2 * phi_02 ** 2 + gamma_x * phi_01 * phi_02
                F_01_ij = integrate.simpson(- 2 * mu * i1 + i2a + i2b + i3, x_grid)
                F_02_ij = integrate.simpson(- 2 * mu * i1 - i2a - i2b - i3, x_grid)
                I1_ij = integrate.simpson(i1, x_grid)
                F_01_i.append(F_01_ij)
                F_02_i.append(F_02_ij)
                I1_i.append(I1_ij)
            F_01.append(F_01_i)
            F_02.append(F_02_i)
            I1.append(I1_i)
        DX = sparse_D(NX, dX)
        DY = sparse_D(NY, dY)

        DY_F1 = np.transpose(Der(DY, np.transpose(F_01)))
        DX_F2 = Der(DX, F_02)
        H = 2 * Der(DX, np.transpose(Der(DY, np.transpose(I1))))
        G1 = np.array(DY_F1 / H)
        G2 = np.array(DX_F2 / H)

        G1_10 = Der(DX, G1)
        G1_01 = np.transpose(Der(DY, np.transpose(G1)))
        G2_10 = Der(DX, G2)
        G2_01 = np.transpose(Der(DY, np.transpose(G2)))

        G1_20 = Der(DX, G1_10)
        G1_11 = Der(DX, G1_01)
        G1_02 = np.transpose(Der(DY, np.transpose(G1_01)))
        G2_20 = Der(DX, G2_10)
        G2_11 = Der(DX, G2_01)
        G2_02 = np.transpose(Der(DY, np.transpose(G2_01)))

        G1_30 = Der(DX, G1_20)
        G1_21 = Der(DX, G1_11)
        G1_12 = np.transpose(Der(DY, np.transpose(G1_11)))
        G1_03 = np.transpose(Der(DY, np.transpose(G1_02)))
        G2_30 = Der(DX, G2_20)
        G2_21 = Der(DX, G2_11)
        G2_12 = np.transpose(Der(DY, np.transpose(G2_11)))
        G2_03 = np.transpose(Der(DY, np.transpose(G2_02)))

        D1_00 = G1[I, J]
        D2_00 = G2[I, J]

        D1_10 = G1_10[I, J]
        D1_01 = G1_01[I, J]
        D2_10 = G2_10[I, J]
        D2_01 = G2_01[I, J]

        M = [[D1_10, D1_01], [D2_10, D2_01]]
        eigenvalues, eigenvectors = np.linalg.eig(M)
        for i in range(len(eigenvalues)):
            plt.scatter(nu, np.real(eigenvalues[i]), c="b")
            plt.scatter(nu, np.imag(eigenvalues[i]), c="r")

        D1_20 = G1_20[I, J]
        D1_11 = G1_11[I, J]
        D1_02 = G1_02[I, J]
        D2_20 = G2_20[I, J]
        D2_11 = G2_11[I, J]
        D2_02 = G2_02[I, J]

        D1_30 = G1_30[I, J]
        D1_21 = G1_21[I, J]
        D1_12 = G1_12[I, J]
        D1_03 = G1_03[I, J]
        D2_30 = G2_30[I, J]
        D2_21 = G2_21[I, J]
        D2_12 = G2_12[I, J]
        D2_03 = G2_03[I, J]

        D1 = [D1_00, D1_10, D1_01, D1_20, D1_11, D1_02, D1_30, D1_21, D1_12, D1_03]
        D2 = [D2_00, D2_10, D2_01, D2_20, D2_11, D2_02, D2_30, D2_21, D2_12, D2_03]
        parameters = np.array([alpha, beta, mu, nu, sigma, gamma])

        alpha_str = f"{alpha:.{3}f}"
        beta_str = f"{beta:.{3}f}"
        mu_str = f"{mu:.{3}f}"
        nu_str = f"{nu:.{3}f}"
        sigma_str = f"{sigma:.{2}f}"
        gamma_str = f"{gamma:.{3}f}"
        savefile = 'D:/mnustes_science/simulation_data/FD/soliton_reduced/gaussian_test/alpha=' + alpha_str + '/beta=' + beta_str + '/mu=' + mu_str + '/nu=' + nu_str + '/gamma=' + gamma_str + '/sigma=' + sigma_str
        if not os.path.exists(savefile):
            os.makedirs(savefile)

        np.savetxt(savefile + '/D1.txt', D1, delimiter=',')
        np.savetxt(savefile + '/D2.txt', D2, delimiter=',')
        np.savetxt(savefile + '/parameters.txt', parameters, delimiter=',')
    plt.show()









