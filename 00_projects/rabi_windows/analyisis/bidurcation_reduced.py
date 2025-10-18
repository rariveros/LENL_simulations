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

    [alpha, beta, mu, nu, gamma_0] = params
    for i in range(len(distances)):
        d = distances[i]
        print("Dist = " + str(d))

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
        Nt_initial = int(0.8 * Nt)
        T = T[Nt_initial:] - T[Nt_initial]

        U1 = UR_Rs[i, Nt_initial:] + 1j * UR_Is[i, Nt_initial:]
        U2 = UL_Rs[i, Nt_initial:] + 1j * UL_Is[i, Nt_initial:]

        R = np.mean(np.abs(U1))
        P = np.mean(np.abs(U2))
        #[R, P] = np.sort(np.array([R, P]))

        phi1 = np.arange(-1.1 * np.pi, 1.1 * np.pi, 0.005)
        phi2 = np.arange(-1.1 * np.pi, 1.1 * np.pi, 0.005)
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

        minima_coords = peak_local_max(-np.sqrt(F1 ** 2 + F2 ** 2), min_distance=10)
        v = np.sqrt(F1 ** 2 + F2 ** 2)[minima_coords[:, 0], minima_coords[:, 1]]
        dot = np.dot(v, v)
        if minima_coords.size != 0 and dot < 0.0001:
            plt.scatter([d] * len(phi1[minima_coords[:, 0]]), phi1[minima_coords[:, 0]], c="b", s=10, edgecolor="k", zorder=10)
            plt.scatter([d] * len(phi2[minima_coords[:, 1]]), phi2[minima_coords[:, 1]], c="r", s=10, edgecolor="k", zorder=10)
    plt.xlabel("$d$", fontsize=15)
    plt.ylabel(r"$\phi_{*}, \theta_{*}$", fontsize=15)
    plt.ylim(-np.pi, np.pi)
    plt.savefig("bifurcation_reduced_flip.png", dpi=300)


