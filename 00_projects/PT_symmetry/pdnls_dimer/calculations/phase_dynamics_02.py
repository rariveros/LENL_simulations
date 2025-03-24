import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

if __name__ == '__main__':
    # Define the grid for the streamplot
    phi1 = np.linspace(-np.pi, np.pi, 250)
    phi2 = np.linspace(-np.pi, np.pi, 250)
    PHI1, PHI2 = np.meshgrid(phi1, phi2)


    # Define the functions F and G
    def F(phi1, phi2, gamma, nu, mu, k, R, P):
        return - nu - k * (P / R) * np.cos(phi1 - phi2) - R ** 2 - gamma * np.sin(2 * phi1)


    def G(phi1, phi2, gamma, nu, mu, k, R, P):
        return - nu - k * (R / P) * np.cos(phi1 - phi2) - P ** 2 + gamma * np.sin(2 * phi2)


    # Fixed Parameters
    nu, mu, k = 0.0, 0.1, 0.05


    # Function to compute R and P based on gamma
    def compute_R_P(gamma):
        R0 = P0 = np.sqrt(-nu + np.sqrt(gamma ** 2 - mu ** 2))
        R1 = (R0 / gamma) * ((np.sqrt(gamma ** 2 - mu ** 2) + mu) / (- 2 * nu + 4 * np.sqrt(gamma ** 2 - mu ** 2)))
        P1 = -(P0 / gamma) * ((np.sqrt(gamma ** 2 - mu ** 2) + mu) / (- 2 * nu + 4 * np.sqrt(gamma ** 2 - mu ** 2)))
        R = R0 + 0 * k * R1  # np.sqrt(-nu + np.sqrt(OMEGAS ** 2 - (mu - k) ** 2))
        P = P0 + 0 * k * P1  # np.sqrt(-nu + np.sqrt(OMEGAS ** 2 - (mu + k) ** 2))
        return R, P


    # Compute values for gamma = 0.15 (Left plot)
    gamma1 = 0.15
    R1, P1 = compute_R_P(gamma1)
    F_values1 = F(PHI1, PHI2, gamma1, nu, mu, k, R1, P1)
    G_values1 = G(PHI1, PHI2, gamma1, nu, mu, k, R1, P1)

    # Compute values for gamma = 0.35 (Right plot)
    gamma2 = 0.25
    R2, P2 = compute_R_P(gamma2)
    F_values2 = F(PHI1, PHI2, gamma2, nu, mu, k, R2, P2)
    G_values2 = G(PHI1, PHI2, gamma2, nu, mu, k, R2, P2)

    # Create figure with 2 subplots
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))  # 1 row, 2 columns

    # Left subplot: Streamplot for gamma = 0.15
    axes[0].streamplot(PHI1, PHI2, F_values1, G_values1, color='k', density=2, linewidth=1)
    contour_G = axes[0].contourf(PHI1, PHI2, np.sqrt(F_values1 ** 2 + G_values1 ** 2), levels=100, cmap='coolwarm')
    fig.colorbar(contour_G, ax=axes[0], label=r"$G(\phi_1, \phi_2)$")
    axes[0].set_xlabel(r'$\phi$')
    axes[0].set_ylabel(r'$\theta$')
    axes[0].set_title(r'Streamplot $(F, G)$ for $\gamma = 0.15$')

    # Right subplot: Streamplot for gamma = 0.35
    axes[1].streamplot(PHI1, PHI2, F_values2, G_values2, color='k', density=2, linewidth=1)
    contour_G = axes[1].contourf(PHI1, PHI2, np.sqrt(F_values2 ** 2 + G_values2 ** 2), levels=100, cmap='coolwarm')
    fig.colorbar(contour_G, ax=axes[1], label=r"$G(\phi_1, \phi_2)$")
    axes[1].set_xlabel(r'$\phi$')
    axes[1].set_ylabel(r'$\theta$')
    axes[1].set_title(r'Streamplot $(F, G)$ for $\gamma = 0.35$')

    # Adjust layout and display
    plt.tight_layout()
    plt.show()
