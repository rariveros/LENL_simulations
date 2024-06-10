from functions import *
from back_process import *
from time_integrators import *
from scipy.signal import argrelextrema

if __name__ == '__main__':
    x_grid = np.arange(-100, 100, 0.1)
    dx = x_grid[1] - x_grid[0]
    alpha = 26.096
    beta = 1
    sigma_i = 20
    nus = [-0.05]#np.arange(0.05, 0.11, 0.005)
    mu = 0.075
    gammas = np.arange(0.15, 0.23, 0.005)

    X = np.arange(-10, 10, 0.1)
    Y = np.arange(-10, 10, 0.1)
    Xs = []
    Ys = []
    nus_scat = []
    for gamma_0 in gammas:
        nu = nus[0]
        print(gamma_0)
        F = []
        G = []
        H = []
        for i in range(len(X)):
            Fi = []
            Gi = []
            Hi = []
            for j in range(len(Y)):
                gammaR_ij = gamma_0 * (np.exp(-X[i] ** 2 / (2 * sigma_i ** 2)))
                gammaI_ij = gamma_0 * (np.exp(-Y[j] ** 2 / (2 * sigma_i ** 2)))
                phaseR_ij = 0.5 * np.arccos(mu / gammaR_ij)
                phaseI_ij = 0.5 * np.arccos(mu / gammaI_ij)
                deltaR_ij = (- nu + np.sqrt(gammaR_ij ** 2 - mu ** 2))
                deltaI_ij = (- nu + np.sqrt(gammaI_ij ** 2 - mu ** 2))
                AR_i = np.sqrt(2 * deltaR_ij) * np.cos(phaseR_ij)
                AI_i = - np.sqrt(2 * deltaI_ij) * np.sin(phaseI_ij)
                phi_01_i = (AR_i / np.cosh(np.sqrt(deltaR_ij / alpha) * (x_grid - X[i])))
                phi_02_i = (AI_i / np.cosh(np.sqrt(deltaI_ij / alpha) * (x_grid - Y[j])))
                Dphi_01_i = np.append(np.diff(phi_01_i) / dx, 0)
                Dphi_02_i = np.append(np.diff(phi_02_i) / dx, 0)
                f_ij = mu * phi_01_i * phi_02_i \
                       - (- (alpha / 2) * Dphi_02_i ** 2 \
                       + (beta / 4) * phi_02_i ** 4 \
                       + (nu / 2) * phi_02_i ** 2\
                       + gamma_0 * (np.exp(-x_grid ** 2 / (2 * sigma_i ** 2))) * phi_01_i * phi_02_i
                          + (beta / 2) * phi_01_i ** 2 * phi_02_i ** 2)
                F_ij = -integrate.simpson(f_ij, x_grid)
                g_ij =  mu * phi_01_i * phi_02_i \
                       + (- (alpha / 2) * Dphi_01_i ** 2 \
                       + (beta / 4) * phi_01_i ** 4 \
                       + (nu / 2) * phi_01_i ** 2\
                       + gamma_0 * (np.exp(-x_grid ** 2 / (2 * sigma_i ** 2))) * phi_01_i * phi_02_i
                       + (beta / 2) * phi_01_i ** 2 * phi_02_i ** 2)
                G_ij = -integrate.simpson(g_ij, x_grid)
                h_ij = phi_01_i * phi_02_i
                H_ij = integrate.simpson(h_ij)
                Fi.append(F_ij)
                Gi.append(G_ij)
                Hi.append(H_ij)
            Fi = np.array(Fi)
            Gi = np.array(Gi)
            Hi = np.array(Hi)
            F.append(Fi)
            G.append(Gi)
            H.append(Hi)
        F = np.array(F)
        G = np.array(G)
        H = np.array(H)
        '''
        plt.figure(figsize=(5, 5))
        plt.pcolormesh(Y, X, F, shading='auto', cmap='jet')
        plt.colorbar(label='Heatmap')
        plt.gca().invert_yaxis()
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.show()
        plt.close()

        plt.figure(figsize=(5, 5))
        plt.pcolormesh(Y, X, G, shading='auto', cmap='jet')
        plt.colorbar(label='Heatmap')
        plt.gca().invert_yaxis()
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.show()
        plt.close()
        '''
        NX = len(X)
        dX = X[1] - X[0]
        NY = len(Y)
        dY = Y[1] - Y[0]
        DX = sparse_D(NX, dX)
        DY = sparse_D(NY, dY)
        dF = np.transpose(Der(DY, np.transpose(F)))
        dG = Der(DX, G)
        ddH = Der(DX, np.transpose(Der(DY, np.transpose(H))))

        V1 = dF[1:-1, 1:-1] #/ ddH[1:-1, 1:-1]
        V2 = dG[1:-1, 1:-1] #/ ddH[1:-1, 1:-1]
        mean_zero_01 = np.mean(np.diff(V1))
        mean_zero_02 = np.mean(np.diff(V2))
        X = X[1:-1]
        Y = Y[1:-1]
        mins_j = []
        for i in range(len(X)):
            for j in range(len(Y)):
                if np.abs(V1[i, j]) < 0.00005 and np.abs(V2[i, j]) < 0.00005:
                    Xs.append(X[i])
                    Ys.append(Y[j])
                    nus_scat.append(gamma_0)

    plt.scatter(nus_scat, Xs, color="b")
    plt.scatter(nus_scat, Ys, color="r")
    plt.grid(alpha=0.2)
    plt.xlabel("$\gamma$", fontsize=20)
    plt.ylabel("$X_s$", fontsize=20)
    plt.show()