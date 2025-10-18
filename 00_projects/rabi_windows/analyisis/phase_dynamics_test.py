import matplotlib.pyplot as plt
from skimage.feature import peak_local_max
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
    d = 17.6
    i = np.argmin(np.abs(distances - d))
    print(i)

    [alpha, beta, mu, nu, gamma_0] = params

    print(np.real(coeffs[i, :]))
    print(np.imag(coeffs[i, :]))

    S11 = coeffs[i, 0]
    S12 = coeffs[i, 1]
    S21 = coeffs[i, 2]
    S22 = coeffs[i, 3]

    P11 = coeffs[i, 4]
    P12 = coeffs[i, 5]
    P21 = coeffs[i, 6]
    P22 = coeffs[i, 7]

    D11 = coeffs[i, 8]
    D21 = coeffs[i, 9]
    D31 = coeffs[i, 10]
    D41 = coeffs[i, 11]
    D51 = coeffs[i, 12]
    D61 = coeffs[i, 13]

    D12 = coeffs[i, 14]
    D22 = coeffs[i, 15]
    D32 = coeffs[i, 16]
    D42 = coeffs[i, 17]
    D52 = coeffs[i, 18]
    D62 = coeffs[i, 19]

    Nt = len(T)
    Nt_initial = int(0.2 * Nt)
    T = T[Nt_initial:] - T[Nt_initial]

    U1 = UR_Rs[i, Nt_initial:] + 1j * UR_Is[i, Nt_initial:]
    U2 = UL_Rs[i, Nt_initial:] + 1j * UL_Is[i, Nt_initial:]

    R = np.mean(np.abs(U1))
    P = np.mean(np.abs(U2))

    phi1 = np.arange(-np.pi, np.pi, 0.005)
    phi2 = np.arange(-np.pi, np.pi, 0.005)
    PHI1, PHI2 = np.meshgrid(phi1, phi2)
    Z1 = R * np.exp(1j * PHI1)
    Z2 = P * np.exp(1j * PHI2) #Hay algo raro con las fases, comparar con la medida y revisar ansatz
    Z1_conj = np.conjugate(Z1)
    Z2_conj = np.conjugate(Z2)
    Z1_mod = np.abs(Z1)
    Z2_mod = np.abs(Z2)

    f1 = - (mu + 1j * nu) * Z1 - 1j * (S11 * Z1 + S12 * Z2) + P11 * Z1_conj + P12 * Z2_conj - 1j * beta * \
        (D11 * (Z2_mod ** 2) * Z2
         + D21 * (Z2_mod ** 2) * Z1
         + D31 * (Z1_mod ** 2) * Z2
         + D41 * (Z1_mod ** 2) * Z1
         + D51 * (Z2 ** 2) * Z1_conj
         + D61 * (Z1 ** 2) * Z2_conj
         )
    f2 = - (mu + 1j * nu) * Z2 - 1j * (S21 * Z1 + S22 * Z2) + P21 * Z1_conj + P22 * Z2_conj - 1j * beta * \
        (D12 * (Z2_mod ** 2) * Z2
         + D22 * (Z2_mod ** 2) * Z1
         + D32 * (Z1_mod ** 2) * Z2
         + D42 * (Z1_mod ** 2) * Z1
         + D52 * (Z2 ** 2) * Z1_conj
         + D62 * (Z1 ** 2) * Z2_conj
         )
    F1 = (1 / R) * np.imag(f1 * np.exp(-1j * PHI1))
    F2 = (1 / P) * np.imag(f2 * np.exp(-1j * PHI2))

    #F1 = np.real((1 / R) * (a10 * R + a01 * P + a30 * R ** 3 + a21 * R ** 2 * P + a12 * R * P ** 2 + a03 * P ** 3))
    #F2 = np.real((1 / P) * (b10 * R + b01 * P + b30 * R ** 3 + b21 * R ** 2 * P + b12 * R * P ** 2 + b03 * P ** 3))

    # Create figure with 2 subplots
    fig, axes = plt.subplots(1, 1, figsize=(4, 4))  # 1 row, 2 columns
    ticks = [0, 0.1, 0.2, 0.3, 0.4]
    #ola

    # Left subplot: Streamplot for gamma = 0.15
    axes.streamplot(PHI1, PHI2, F1, F2, color='k', density=2.0, linewidth=1, arrowsize=1)
    contour_G1 = axes.contourf(PHI1, PHI2, np.sqrt(F1 ** 2 + F2 ** 2), levels=100, cmap='turbo', alpha=1.0) #np.sqrt(F1 ** 2 + F2 ** 2)
    minima_coords = peak_local_max(-np.log(np.sqrt(F1 ** 2 + F2 ** 2)), min_distance=5)
    print(minima_coords)
    plt.scatter(phi1[minima_coords[:, 1]], phi2[minima_coords[:, 0]], c="w", s=25, edgecolor="k", zorder=10)

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
