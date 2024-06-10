from functions import *
from back_process import *
from time_integrators import *
from scipy.signal import argrelextrema
from skimage import measure
from shapely.geometry import LineString

if __name__ == '__main__':
    x_grid = np.arange(-75, 75, 0.1)
    dx = x_grid[1] - x_grid[0]
    alpha = 26.096
    beta = 1
    sigma_i = 19
    nu = -0.1
    mu = 0.075
    gamma_0 = 0.22

    x = X = np.arange(-20, 20, 0.1)
    y = Y = np.arange(-20, 20, 0.1)
    NX = len(X)
    dX = X[1] - X[0]
    NY = len(Y)
    dY = Y[1] - Y[0]
    DX = sparse_D(NX, dX)
    DY = sparse_D(NY, dY)
    X, Y = np.meshgrid(X, Y)
    gammaR = gamma_0 * (np.exp(-X ** 2 / (2 * sigma_i ** 2)))
    gammaI = gamma_0 * (np.exp(-Y ** 2 / (2 * sigma_i ** 2)))
    phaseR = 0.5 * np.arccos(mu / gammaR)
    phaseI = 0.5 * np.arccos(mu / gammaI)
    deltaR = (- nu + np.sqrt(gammaR ** 2 - mu ** 2))
    deltaI = (- nu + np.sqrt(gammaI ** 2 - mu ** 2))
    AR = np.sqrt(2 * deltaR) * np.cos(phaseR)
    AI = - np.sqrt(2 * deltaI) * np.sin(phaseI)
    bR = np.sqrt(deltaR / alpha)
    bI = np.sqrt(deltaI / alpha)
    dX = Y - X
    A = 0.5 * (AR * AI / bR) * (2 - (bR * bI * dX ** 2 / 3))
    B1 = -alpha * AR ** 2 * bR / 3 + nu * AR ** 2 / bR + 2 * beta * AR ** 4 / (3 * bR)
    B2 = -alpha * AI ** 2 * bI / 3 + nu * AI ** 2 / bI + 2 * beta * AI ** 4 / (3 * bI)
    C = (beta * AR * AI / (2 * bR)) * (4/3 + bR ** 2 * dX ** 2 * np.pi / 8) + (gamma_0 * AR * AI / bR) * np.exp(-Y * X / (2 * sigma_i ** 2)) * (2 - bR ** 2 * dX ** 2 / 3 - X * dX / (sigma_i ** 2))

    F = - 2 * mu * A + C #+ B1
    G = - 2 * mu * A - C #- B2

    dF = Der(DY, F)
    dG = np.transpose(Der(DX, np.transpose(G)))
    V1 = dF[1:-1, 1:-1]  # / ddH[1:-1, 1:-1]
    V2 = dG[1:-1, 1:-1]
    x = x[1:-1]
    y = y[1:-1]

    fig, (ax1, ax2) = plt.subplots(1, 2)
    ax1.pcolormesh(x, y, np.transpose(V1), shading='auto', cmap='jet')
    ax1.set_label("X")
    ax1.set_ylabel("Y")
    ax1.contour(x, y, np.transpose(V1), [0])

    ax2.pcolormesh(x, y, np.transpose(V2), shading='auto', cmap='jet')
    ax2.set_xlabel("X")
    ax2.contour(x, y, np.transpose(V2), [0])

    plt.show()
    plt.close()
