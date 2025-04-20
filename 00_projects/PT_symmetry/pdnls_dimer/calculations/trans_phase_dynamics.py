import matplotlib.pyplot as plt
from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    # Define the grid
    phi1 = np.arange(-np.pi, np.pi, 0.005)
    phi2 = np.arange(-np.pi, np.pi, 0.005)
    PHI1, PHI2 = np.meshgrid(phi1, phi2)

    # Define the functions F and G
    def F(phi1, phi2, gamma, nu, mu, k):
        R = np.sqrt(((gamma + mu) / gamma) * (-nu + np.sqrt(gamma ** 2 - mu ** 2)))
        P = np.sqrt(R ** 2 * (gamma - mu) / (gamma + mu))
        delta_1 = R ** 2 * np.cos(phi1) ** 2 + P ** 2 * np.sin(phi2) ** 2
        delta_2 = R ** 2 * np.sin(phi1) ** 2 + P ** 2 * np.cos(phi2) ** 2
        return - k - nu * (P / R) * np.cos(phi1 - phi2) - (P / R) * np.cos(phi1) * np.cos(phi2) * delta_2 - (P / R) * np.sin(phi1) * np.sin(phi2) * delta_1


    def G(phi1, phi2, gamma, nu, mu, k):
        R = np.sqrt(((gamma + mu) / gamma) * (-nu + np.sqrt(gamma ** 2 - mu ** 2)))
        P = np.sqrt(R ** 2 * (gamma - mu) / (gamma + mu))
        delta_1 = R ** 2 * np.cos(phi1) ** 2 + P ** 2 * np.sin(phi2) ** 2
        delta_2 = R ** 2 * np.sin(phi1) ** 2 + P ** 2 * np.cos(phi2) ** 2
        return - k - nu * (R / P) * np.cos(phi1 - phi2) - (R / P) * np.cos(phi1) * np.cos(phi2) * delta_1 - (R / P) * np.sin(phi1) * np.sin(phi2) * delta_2


    # Parameters
    [nu, mu, k] = [0.1, 0.1, 0.05]

    gamma1 = 0.195
    F_values1 = F(PHI1, PHI2, gamma1, nu, mu, k)
    G_values1 = G(PHI1, PHI2, gamma1, nu, mu, k)

    # Compute values for gamma = 0.35 (Right plot)
    gamma2 = 0.25
    F_values2 = F(PHI1, PHI2, gamma2, nu, mu, k)
    G_values2 = G(PHI1, PHI2, gamma2, nu, mu, k)

    # Create figure with 2 subplots
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))  # 1 row, 2 columns

    # Left subplot: Streamplot for gamma = 0.15
    axes[0].streamplot(PHI1, PHI2, F_values1, G_values1, color='k', density=4, linewidth=1)
    contour_G = axes[0].contourf(PHI1, PHI2, np.sqrt(F_values1 ** 2 + G_values1 ** 2), levels=100, cmap='turbo', alpha=1.0)
    cbar = fig.colorbar(contour_G, ax=axes[0])
    cbar.ax.tick_params(labelsize=15)
    cbar.set_label(r"$|\vec{F}(\phi_1, \phi_2)|$", rotation=0, size=18, labelpad=-50, y=1.1)
    axes[0].set_xlabel(r'$\phi$', fontsize=20)
    axes[0].set_ylabel(r'$\theta$', fontsize=20)
    axes[0].set_title(r'$\gamma = 0.195$', fontsize=15)
    axes[0].tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
    axes[0].tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=True)
    #axes[0].set_xlim(-np.pi, np.pi)

    # Right subplot: Streamplot for gamma = 0.35
    axes[1].streamplot(PHI1, PHI2, F_values2, G_values2, color='k', density=4, linewidth=1)
    contour_G = axes[1].contourf(PHI1, PHI2, np.sqrt(F_values2 ** 2 + G_values2 ** 2), levels=100, cmap='turbo', alpha=1.0)
    cbar = fig.colorbar(contour_G, ax=axes[1])
    cbar.ax.tick_params(labelsize=15)
    cbar.set_label(r"$|\vec{F}(\phi_1, \phi_2)|$", rotation=0, size=18, labelpad=-30, y=1.1)
    axes[1].set_xlabel(r'$\phi$', fontsize=20)
    axes[1].set_ylabel(r'$\theta$', fontsize=20)
    axes[1].set_title(r'$\gamma = 0.250$', fontsize=15)
    axes[1].tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=False, labelright=False)
    axes[1].tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=True)

    # Adjust layout and display
    plt.savefig('trans_phase_dynamics.png', dpi=300)
    plt.tight_layout()
    plt.show()