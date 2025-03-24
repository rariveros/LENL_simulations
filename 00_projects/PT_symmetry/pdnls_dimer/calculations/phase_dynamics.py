import matplotlib.pyplot as plt
from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    # Define the grid
    phi1 = np.arange(-np.pi, np.pi, 0.02)
    phi2 = np.arange(-np.pi, np.pi, 0.02)
    PHI1, PHI2 = np.meshgrid(phi1, phi2)


    # Define the functions F and G
    def F(phi1, phi2, gamma, nu, mu, k, R, P):
        return - nu - k * (P / R) * np.cos(phi1 - phi2) - R ** 2 - gamma * np.sin(2 * phi1)


    def G(phi1, phi2, gamma, nu, mu, k, R, P):
        return - nu - k * (R / P) * np.cos(phi1 - phi2) - P ** 2 + gamma * np.sin(2 * phi2)


    # Parameters
    [gamma, nu, mu, k] = [0.15, 0.1, 0.1, 0.05]
    R0 = P0 = np.sqrt(-nu + np.sqrt(gamma ** 2 - mu ** 2))
    R1 = (R0 / gamma) * ((np.sqrt(gamma ** 2 - mu ** 2) + mu) / (- 2 * nu + 4 * np.sqrt(gamma ** 2 - mu ** 2)))
    P1 = -(P0 / gamma) * ((np.sqrt(gamma ** 2 - mu ** 2) + mu) / (- 2 * nu + 4 * np.sqrt(gamma ** 2 - mu ** 2)))
    R = R0 + k * R1
    P = P0 + k * P1

    # Compute function values
    F_values = F(PHI1, PHI2, gamma, nu, mu, k, R, P)
    G_values = G(PHI1, PHI2, gamma, nu, mu, k, R, P)

    # Create figure with 2 subplots
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))  # 1 row, 2 columns

    # Plot F(phi1, phi2)
    contour_F = axes[0].contourf(PHI1, PHI2, F_values, levels=50, cmap='coolwarm')
    fig.colorbar(contour_F, ax=axes[0], label=r"$F(\phi_1, \phi_2)$")
    axes[0].contour(PHI1, PHI2, F_values, levels=[0], colors='black', linewidths=2)
    axes[0].contour(PHI1, PHI2, G_values, levels=[0], colors='gray', linewidths=2)
    axes[0].set_xlabel(r'$\phi_{-}$')
    axes[0].set_ylabel(r'$\phi_{+}$')
    axes[0].set_title(r'$F(\phi_1, \phi_2)$')

    # Plot G(phi1, phi2)
    contour_G = axes[1].contourf(PHI1, PHI2, G_values, levels=50, cmap='coolwarm')
    fig.colorbar(contour_G, ax=axes[1], label=r"$G(\phi_1, \phi_2)$")
    axes[1].contour(PHI1, PHI2, F_values, levels=[0], colors='gray', linewidths=2)
    axes[1].contour(PHI1, PHI2, G_values, levels=[0], colors='black', linewidths=2)
    axes[1].set_xlabel(r'$\phi_{-}$')
    axes[1].set_ylabel(r'$\phi_{+}$')
    axes[1].set_title(r'$G(\phi_1, \phi_2)$')

    # Adjust layout and display
    plt.tight_layout()
    plt.show()