import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

if __name__ == '__main__':
    # === parámetros ===
    a = 1.0/6 # ancho del solitón
    sigma = 15.0  # ancho de la gaussiana
    c = 1.0  # fase base
    d = 0.0  # frecuencia espacial
    b_vals = np.arange(-50, 50, 0.1)  # desplazamiento del solitón


    # === integrando ===
    def integrand(x, b, a, sigma, c, d):
        return np.exp(-x ** 2 / (2 * sigma ** 2)) * (1 / np.cosh(a * (x - b))) ** 2 * np.sin(2 * (c - d * (x - b)))


    # === integral numérica para cada b ===
    I_vals = []
    for b in b_vals:
        val, _ = quad(integrand, -np.inf, np.inf, args=(b, a, sigma, c, d))
        I_vals.append(val)

    # === graficar ===
    plt.figure(figsize=(7, 4))
    plt.plot(b_vals, I_vals, color='r', lw=2)
    plt.plot(b_vals, np.exp(-b_vals ** 2 / (2 * sigma ** 2)), color='g', lw=2)
    plt.plot(b_vals, (1 / np.cosh(a * b_vals)), color='b', lw=2)

    plt.xlabel("b (desplazamiento)", fontsize=12)
    plt.ylabel("Convolución I(b)", fontsize=12)
    plt.title(r"$I(b)=\int e^{-x^2/(2\sigma^2)}\,\mathrm{sech}^2[a(x-b)]\,\sin[2(c-d(x-b))]\,dx$", fontsize=11)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()