from functions import *
from back_process import *
from time_integrators import *

def phis(alpha, beta, nu, mu, gamma, sigma, X, Y, x_grid, dx):
    gamma_01 = gamma * np.exp(- X ** 2 / (2 * sigma ** 2))
    gamma_02 = gamma * np.exp(- Y ** 2 / (2 * sigma ** 2))
    delta_01 = - nu + np.sqrt(gamma_01 ** 2 - mu ** 2)
    delta_02 = - nu + np.sqrt(gamma_02 ** 2 - mu ** 2)
    theta_01 = 0.5 * np.arccos(mu / gamma_01)
    theta_02 = 0.5 * np.arccos(mu / gamma_02)
    phi_01 = (1 / (np.cosh(np.sqrt(delta_01 / alpha) * (x_grid - X)))) * np.sqrt(2 * delta_01) * np.cos(theta_01) * (beta ** (-2))
    phi_02 = - (1 / (np.cosh(np.sqrt(delta_02 / alpha) * (x_grid - Y)))) * np.sqrt(2 * delta_02) * np.sin(theta_02) * (beta ** (-2))
    Dphi_01 = np.append(np.diff(phi_01) / dx, 0)
    Dphi_02 = np.append(np.diff(phi_02) / dx, 0)
    return [phi_01, phi_02, Dphi_01, Dphi_02]


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
    EIGS = []
    for nu in nus:
        print("nu = {}".format(nu))
        [alpha, beta, nu, mu, gamma, sigma] = [6.524, 1, nu, 0.075, 0.18, 15]
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








