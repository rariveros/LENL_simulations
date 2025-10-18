import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    frequencies = []
    powers = []
    modules = []
    disc = 'D:/'

    initial_dir_data = str(disc) + 'Users/mnustes_science/PT_fluids/mnustes_science/simulation_data'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    distances = np.loadtxt(directory + '/analysis/dists.txt', delimiter=',')
    T = np.loadtxt(directory + '/analysis/t_grid.txt', delimiter=',')
    UR_Rs = np.loadtxt(directory + '/analysis/U1_Rs.txt', delimiter=',')
    UR_Is = np.loadtxt(directory + '/analysis/U1_Is.txt', delimiter=',')
    UL_Rs = np.loadtxt(directory + '/analysis/U2_Rs.txt', delimiter=',')
    UL_Is = np.loadtxt(directory + '/analysis/U2_Is.txt', delimiter=',')
    coeffs = np.loadtxt(directory + '/analysis/coefs.txt', delimiter=',', dtype=complex)
    params = np.loadtxt(directory + '/analysis/params.txt', delimiter=',', dtype=complex)

    save_directory = directory + '/analysis'
    d = 40
    i = np.argmin(np.abs(distances - d))
    print(i)

    [alpha, beta, mu, nu, gamma_0] = params

    print(np.real(coeffs[i, :]))
    print(np.imag(coeffs[i, :]))

    Sigma_11 = coeffs[i, 0]
    Sigma_12 = coeffs[i, 1]
    Sigma_21 = coeffs[i, 2]
    Sigma_22 = coeffs[i, 3]

    Pi_11 = coeffs[i, 4]
    Pi_12 = coeffs[i, 5]
    Pi_21 = coeffs[i, 6]
    Pi_22 = coeffs[i, 7]

    Delta_11 = coeffs[i, 8]
    Delta_21 = coeffs[i, 9]
    Delta_31 = coeffs[i, 10]
    Delta_41 = coeffs[i, 11]
    Delta_51 = coeffs[i, 12]
    Delta_61 = coeffs[i, 13]

    Delta_12 = coeffs[i, 14]
    Delta_22 = coeffs[i, 15]
    Delta_32 = coeffs[i, 16]
    Delta_42 = coeffs[i, 17]
    Delta_52 = coeffs[i, 18]
    Delta_62 = coeffs[i, 19]

    Nt = len(T)
    Nt_initial = int(0.2 * Nt)
    T = T[Nt_initial:] - T[Nt_initial]

    UR = UR_Rs[i, Nt_initial:] + 1j * UR_Is[i, Nt_initial:]
    UL = UL_Rs[i, Nt_initial:] + 1j * UL_Is[i, Nt_initial:]
    U1 = 0.5 * (UL + 1j * UR)
    U2 = 0.5 * (UL - 1j * UR)

    R = np.abs(U1)[-1]
    P = np.abs(U2)[-1]

    phi1 = np.arange(-np.pi, np.pi, 0.02)
    phi2 = np.arange(-np.pi, np.pi, 0.02)
    PHI1, PHI2 = np.meshgrid(phi1, phi2)

    a10 = - nu - np.real(Sigma_11)
    a01 = np.cos(PHI1 + PHI2) * np.imag(Pi_12) - np.sin(PHI1 + PHI2) * np.real(Pi_12)
    a30 = 0
    a21 = -np.cos(-(PHI1 - PHI2)) * np.imag(Delta_41) + np.sin(-(PHI1 - PHI2)) * np.real(Delta_41)
    a12 = -np.cos(2 * (PHI1 - PHI2)) * np.imag(Delta_51) - np.sin(2 * (PHI1 - PHI2)) * np.real(Delta_51)
    a03 = -np.cos(-(PHI1 - PHI2)) * np.imag(Delta_21) + np.sin(-(PHI1 - PHI2)) * np.real(Delta_21)

    b10 = np.cos(PHI1 + PHI2) * np.imag(Pi_21) - np.sin(PHI1 + PHI2) * np.real(Pi_21)
    b01 = - nu - np.real(Sigma_22)
    b30 = -np.cos(-(PHI1 - PHI2)) * np.imag(Delta_32) - np.sin(-(PHI1 - PHI2)) * np.real(Delta_32)
    b21 = -np.cos(2 * (PHI1 - PHI2)) * np.imag(Delta_62) + np.sin(2 * (PHI1 - PHI2)) * np.real(Delta_62)
    b12 = -np.cos(-(PHI1-PHI2)) * np.imag(Delta_12) - np.sin(-(PHI1-PHI2)) * np.real(Delta_12)
    b03 = 0

    F1 = np.real((1 / R) * (a10 * R + a01 * P + a30 * R ** 3 + a21 * R ** 2 * P + a12 * R * P ** 2 + a03 * P ** 3))
    F2 = np.real((1 / P) * (b10 * R + b01 * P + b30 * R ** 3 + b21 * R ** 2 * P + b12 * R * P ** 2 + b03 * P ** 3))

    # Create figure with 2 subplots
    fig, axes = plt.subplots(1, 1, figsize=(4, 4))  # 1 row, 2 columns
    ticks = [0, 0.1, 0.2, 0.3, 0.4]
    #ola

    # Left subplot: Streamplot for gamma = 0.15
    #axes.streamplot(PHI1, PHI2, F1, F2, color='k', density=5, linewidth=0.7, arrowsize=0.7)
    contour_G1 = axes.contourf(PHI1, PHI2, F2, levels=100, cmap='turbo', alpha=1.0) #np.sqrt(F1 ** 2 + F2 ** 2)
    #axes.contour(
    #    PHI1, PHI2, F1,
    #    levels=[1e-3],  # umbral pequeño para evitar log(0)
    #    colors='k', linewidths=1.0
    #)
    axes.contour(
        PHI1, PHI2, F1,
        levels=[1e-3],  # umbral pequeño para evitar log(0)
        colors='k', linewidths=1.0
    )
    cbar = fig.colorbar(contour_G1, ax=axes)#, ticks=ticks)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label(r"$|\vec{F}(\theta, \phi)|$", rotation=0, size=14, labelpad=-20, y=1.23)
    axes.set_xlabel(r'$\phi$', fontsize=14)
    axes.set_ylabel(r'$\theta$', fontsize=14)
    #axes[0].set_title(r'$\kappa = 0.03$', fontsize=12)
    xticks = [-np.pi, 0, np.pi]
    xtick_labels = [r"$-\pi$", r"$0$", r"$\pi$"]
    axes.set_xticks(xticks)
    axes.set_xticklabels(xtick_labels)
    axes.set_yticks(xticks)
    axes.set_yticklabels(xtick_labels)
    axes.tick_params(axis="y", direction="in", labelsize=12, left=True, right=True, labelleft=True, labelright=False)
    axes.tick_params(axis="x", direction="in", labelsize=12, top=True, bottom=True, labeltop=False, labelbottom=True)
    axes.set_xlim(-np.pi, np.pi)
    axes.set_ylim(-np.pi, np.pi)

    fig.subplots_adjust(left=0.1, right=0.9, bottom=0.25, top=0.85)
    plt.savefig('reduced_phase_dynamics.png', dpi=300)
    plt.close()
