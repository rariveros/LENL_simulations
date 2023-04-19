from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    # Number of sample points
    N = 600
    # sample spacing
    T = 1.0 / 800.0
    x = np.linspace(0.0, N*T, N, endpoint=False)
    f = 10
    y = np.sin(2 * np.pi * f * x)
    yf = fft(y)
    xf = fftfreq(N, T)[:N//2]

    plt.plot(xf, 2.0/N * np.abs(yf[0:N//2]))
    plt.grid()
    plt.show()