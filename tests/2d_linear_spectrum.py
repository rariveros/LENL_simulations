from functions import *
from back_process import *
from time_integrators import *
import scipy.sparse as sp
import scipy.sparse.linalg as spla

if __name__ == '__main__':
    import numpy as np
    import scipy.sparse.linalg as spla
    import matplotlib.pyplot as plt
    from scipy.fft import fft2, ifft2, fftfreq

    # Parámetros físicos solicitados
    mu, nu, alpha = 0.1, 0.1, 1.0
    sigma = 5
    gamma_amp = 0.0

    # Dominio y discretización
    L = 20
    x_min, x_max = -L, L
    y_min, y_max = -L, L
    dx, dy = 0.25, 0.25
    Lx, Ly = x_max - x_min, y_max - y_min
    x = np.arange(x_min, x_max, dx)
    y = np.arange(y_min, y_max, dy)
    Nx, Ny = len(x), len(y)
    X, Y = np.meshgrid(x, y, indexing='ij')
    N = Nx * Ny

    # Espacio de Fourier
    kx = fftfreq(Nx, d=Lx / (2 * np.pi * Nx))
    ky = fftfreq(Ny, d=Ly / (2 * np.pi * Ny))
    KX, KY = np.meshgrid(kx, ky, indexing='ij')
    K2 = KX ** 2 + KY ** 2

    # Campo gamma y su FFT
    gamma_prime = gamma_amp * np.exp(-(X ** 2 + Y ** 2) / (2 * sigma ** 2))
    gamma_fft = fft2(gamma_prime)


    # Operador lineal espectral
    def apply_L(v):
        A = v[:N].reshape((Nx, Ny))
        A_star = v[N:].reshape((Nx, Ny))

        A_fft = fft2(A)
        A_lap = ifft2(-1j * alpha * (-K2) * A_fft)
        L11 = - (mu + 1j * nu) * A + A_lap

        A_star_conv = ifft2(gamma_fft * fft2(A_star))
        L12 = A_star_conv

        A_star_fft = fft2(A_star)
        A_star_lap = ifft2(1j * alpha * K2 * A_star_fft)
        L22 = - (mu - 1j * nu) * A_star + A_star_lap

        L21 = ifft2(gamma_fft * fft2(A))

        top = L11 + L12
        bottom = L21 + L22
        return np.concatenate([top.ravel(), bottom.ravel()])


    # LinearOperator para eigs
    Lop = spla.LinearOperator((2 * N, 2 * N), matvec=apply_L, dtype=np.complex128)

    # Autovalor más inestable
    eigvals, eigvecs = spla.eigs(Lop, k=1, which='LR')
    mode = eigvecs[:, 0]
    A_mode = mode[:N].reshape((Nx, Ny))
    A_real = np.real(A_mode)
    A_imag = np.imag(A_mode)

    # Visualización
    fig, axes = plt.subplots(1, 3, figsize=(12, 5))
    axes[0].imshow(A_real, extent=[x_min, x_max, y_min, y_max], origin='lower', cmap='RdBu')
    axes[0].set_title('Parte real de A')
    axes[1].imshow(A_imag, extent=[x_min, x_max, y_min, y_max], origin='lower', cmap='RdBu')
    axes[1].set_title('Parte imaginaria de A')
    axes[2].imshow(gamma_prime, extent=[x_min, x_max, y_min, y_max], origin='lower', cmap='RdBu')
    axes[2].set_title('Gamma')
    plt.tight_layout()
    plt.show()

    print(f"Autovalor dominante: {eigvals[0]}")
