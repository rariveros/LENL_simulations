import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

if __name__ == '__main__':
    # Define the grid for the streamplot
    r = np.linspace(0.01, 0.5, 250)
    p = np.linspace(0.01, 0.5, 250)
    R, P = np.meshgrid(r, p)


    # Define the functions F and G
    def F1(phi1, phi2, gamma, nu, mu, k, R, P):
        return - mu * R - k * P * np.sin(phi1 - phi2) + gamma * R *np.cos(2 * phi1)

    def G1(phi1, phi2, gamma, nu, mu, k, R, P):
        return - nu - k * (P / R) * np.cos(phi1 - phi2) - R ** 2 - gamma * np.sin(2 * phi1)

    def F2(phi1, phi2, gamma, nu, mu, k, R, P):
        return - mu * P + k * R * np.sin(phi1 - phi2) - gamma * P *np.cos(2 * phi2)

    def G2(phi1, phi2, gamma, nu, mu, k, R, P):
        return - nu - k * (R / P) * np.cos(phi1 - phi2) - P ** 2 + gamma * np.sin(2 * phi2)

    # Fixed Parameters
    nu, mu, k = 0.1, 0.1, 0.05
    gamma = 0.5
    PHI1 = 0.5 * np.arccos(mu / gamma)
    PHI2 = 0.5 * np.arccos(-mu / gamma)

    F_values1 = F1(PHI1, PHI2, gamma, nu, mu, k, R, P)
    G_values1 = G1(PHI1, PHI2, gamma, nu, mu, k, R, P)
    F_values2 = F2(PHI1, PHI2, gamma, nu, mu, k, R, P)
    G_values2 = G2(PHI1, PHI2, gamma, nu, mu, k, R, P)

    # Create figure with 2 subplots
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))  # 1 row, 2 columns

    # Left subplot: Streamplot for gamma = 0.15
    axes[0].streamplot(R, P, F_values1, F_values2, color='k', density=2, linewidth=1)
    contour_G = axes[0].contourf(R, P, np.sqrt(F_values1 ** 2 + F_values2 ** 2), levels=100, cmap='coolwarm')
    fig.colorbar(contour_G, ax=axes[0], label=r"$G(\phi_1, \phi_2)$")
    axes[0].set_xlabel(r'$R$')
    axes[0].set_ylabel(r'$P$')

    # Right subplot: Streamplot for gamma = 0.35
    axes[1].streamplot(R, P, G_values1, G_values2, color='k', density=2, linewidth=1)
    contour_G = axes[1].contourf(R, P, np.sqrt(G_values1 ** 2 + G_values2 ** 2), levels=100, cmap='coolwarm')
    fig.colorbar(contour_G, ax=axes[1], label=r"$G(\phi_1, \phi_2)$")
    axes[1].set_xlabel(r'$R$')
    axes[1].set_ylabel(r'$P$')

    # Adjust layout and display
    plt.tight_layout()
    plt.show()
