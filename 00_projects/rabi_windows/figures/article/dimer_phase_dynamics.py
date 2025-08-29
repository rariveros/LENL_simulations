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
    gamma = 0.2
    k1 = 0.03
    F_values1 = F(PHI1, PHI2, gamma, nu, mu, k1)
    G_values1 = G(PHI1, PHI2, gamma, nu, mu, k1)

    # Compute values for gamma = 0.35 (Right plot)
    k2 = 0.05
    F_values2 = F(PHI1, PHI2, gamma, nu, mu, k2)
    G_values2 = G(PHI1, PHI2, gamma, nu, mu, k2)

    # Create figure with 2 subplots
    fig, axes = plt.subplots(1, 2, figsize=(8, 2))  # 1 row, 2 columns
    ticks = [0, 0.1, 0.2, 0.3, 0.4]

    # Left subplot: Streamplot for gamma = 0.15
    axes[0].streamplot(PHI1, PHI2, F_values1, G_values1, color='k', density=2.2, linewidth=0.7, arrowsize=0.7)
    contour_G1 = axes[0].contourf(PHI1, PHI2, np.sqrt(F_values1 ** 2 + G_values1 ** 2), levels=100, cmap='turbo', alpha=1.0)
    cbar = fig.colorbar(contour_G1, ax=axes[0], ticks=ticks)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label(r"$|\vec{F}(\theta, \phi)|$", rotation=0, size=14, labelpad=-20, y=1.23)
    axes[0].set_xlabel(r'$\phi$', fontsize=14)
    axes[0].set_ylabel(r'$\theta$', fontsize=14)
    axes[0].set_title(r'$\kappa = 0.03$', fontsize=12)
    xticks = [-np.pi, 0, np.pi]
    xtick_labels = [r"$-\pi$", r"$0$", r"$\pi$"]
    axes[0].set_xticks(xticks)
    axes[0].set_xticklabels(xtick_labels)
    axes[0].set_yticks(xticks)
    axes[0].set_yticklabels(xtick_labels)
    axes[0].tick_params(axis="y", direction="in", labelsize=12, left=True, right=True, labelleft=True, labelright=False)
    axes[0].tick_params(axis="x", direction="in", labelsize=12, top=True, bottom=True, labeltop=False, labelbottom=True)
    axes[0].set_xlim(-np.pi, np.pi)

    # Right subplot: Streamplot for gamma = 0.35
    axes[1].streamplot(PHI1, PHI2, F_values2, G_values2, color='k', density=2.2, linewidth=0.7, arrowsize=0.7)
    contour_G2 = axes[1].contourf(PHI1, PHI2, np.sqrt(F_values2 ** 2 + G_values2 ** 2), levels=100, cmap='turbo', alpha=1.0)
    cbar = fig.colorbar(contour_G2, ax=axes[1], ticks=ticks)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label(r"$|\vec{F}(\theta, \phi)|$", rotation=0, size=14, labelpad=-20, y=1.23)
    xticks = [-np.pi, 0, np.pi]
    xtick_labels = [r"$-\pi$", r"$0$", r"$\pi$"]
    axes[1].set_xticks(xticks)
    axes[1].set_xticklabels(xtick_labels)
    axes[1].set_yticks(xticks)
    axes[1].set_yticklabels(xtick_labels)
    axes[1].set_xlabel(r'$\phi$', fontsize=14)
    axes[1].set_title(r'$\kappa = 0.05$', fontsize=12)
    axes[1].tick_params(axis="y", direction="in", labelsize=12, left=True, right=True, labelleft=False, labelright=False)
    axes[1].tick_params(axis="x", direction="in", labelsize=12, top=True, bottom=True, labeltop=False, labelbottom=True)

    # Adjust layout and display
    fig.subplots_adjust(left=0.1, right=0.9, bottom=0.25, top=0.85)
    plt.savefig('dimer_phase_dynamics.png', dpi=300)