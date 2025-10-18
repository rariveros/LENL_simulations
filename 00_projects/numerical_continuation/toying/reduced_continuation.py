import numpy as np
import matplotlib.pyplot as plt

# =========================
#   BLOQUE A — SETUP
# =========================

def get_coefs_at(d, distances, coeffs):
    """Interpolación lineal compleja entre coeficientes adyacentes."""
    idx = np.searchsorted(distances, d) - 1
    idx = np.clip(idx, 0, len(distances)-2)
    d1, d2 = distances[idx], distances[idx+1]
    w = (d - d1) / (d2 - d1)
    return coeffs[idx, :] * (1 - w) + coeffs[idx+1, :] * w


def central_diff(A, dx):
    """Derivada centrada con manejo de bordes."""
    dA = np.zeros_like(A)
    dA[1:-1, :] = (A[2:, :] - A[:-2, :]) / (2 * dx)
    dA[0, :] = (A[1, :] - A[0, :]) / dx
    dA[-1, :] = (A[-1, :] - A[-2, :]) / dx
    return dA


def build_F(PHI1, PHI2, coeffs_row, params):
    """Genera F1 y F2 a partir de los coeficientes (complejos) y parámetros."""
    alpha, beta, mu, nu, gamma_0 = params

    (S11, S12, S21, S22,
     P11, P12, P21, P22,
     D11, D21, D31, D41, D51, D61,
     D12, D22, D32, D42, D52, D62) = coeffs_row

    R, P = 1.0, 1.0

    Z1 = R * np.exp(1j * PHI1)
    Z2 = P * np.exp(1j * PHI2)
    Z1c, Z2c = np.conjugate(Z1), np.conjugate(Z2)
    Z1m, Z2m = np.abs(Z1), np.abs(Z2)

    f1 = -(mu + 1j * nu) * Z1 - 1j * (S11 * Z1 + S12 * Z2) + P11 * Z1c + P12 * Z2c - 1j * beta * (D11 * (Z2m**2) * Z2 + D21 * (Z2m**2) * Z1 + D31 * (Z1m**2) * Z2 +
                        D41 * (Z1m**2) * Z1 + D51 * (Z2**2) * Z1c + D61 * (Z1**2) * Z2c)
    f2 = -(mu + 1j * nu) * Z2 - 1j * (S21 * Z1 + S22 * Z2) + P21 * Z1c + P22 * Z2c - 1j * beta * (D12 * (Z2m**2) * Z2 + D22 * (Z2m**2) * Z1 + D32 * (Z1m**2) * Z2 +
                        D42 * (Z1m**2) * Z1 + D52 * (Z2**2) * Z1c + D62 * (Z1**2) * Z2c)
    F1 = (1 / R) * np.imag(f1 * np.exp(-1j * PHI1))
    F2 = (1 / P) * np.imag(f2 * np.exp(-1j * PHI2))
    return F1, F2


# =========================
#   BLOQUE B — CONTINUACIÓN
# =========================

def continuation(directory, d0, ds=0.1, steps=40):
    """Continuación numérica interpolando coeficientes en cada paso (seguimiento local de rama)."""
    # Carga datos
    distances = np.loadtxt(directory + '/analysis/dists.txt', delimiter=',')
    coeffs = np.loadtxt(directory + '/analysis/coefs.txt', delimiter=',', dtype=complex)
    params = np.loadtxt(directory + '/analysis/params.txt', delimiter=',', dtype=complex)
    alpha, beta, mu, nu, gamma_0 = params

    phi1 = np.arange(-np.pi, np.pi, 0.005)
    phi2 = np.arange(-np.pi, np.pi, 0.005)
    PHI1, PHI2 = np.meshgrid(phi1, phi2)
    dphi = phi1[1] - phi1[0]

    minima = []
    prev_idx = None

    d = d0
    for s in range(steps):
        coeffs_row = get_coefs_at(d, distances, coeffs)
        F1, F2 = build_F(PHI1, PHI2, coeffs_row, params)
        normF = np.sqrt(F1**2 + F2**2)

        # Localiza mínimo de |F|
        if s == 0:
            min_idx = np.unravel_index(np.argmin(normF), normF.shape)
            prev_idx = min_idx
        else:
            dist = (PHI1 - PHI1[prev_idx])**2 + (PHI2 - PHI2[prev_idx])**2
            mask = dist > (0.3)**2
            masked = np.ma.array(normF, mask=mask)
            min_idx = np.unravel_index(np.argmin(masked), masked.shape)
            prev_idx = min_idx

        # --- Refinar el mínimo con Newton local ---
        phi1_guess = phi1[min_idx[1]]
        phi2_guess = phi2[min_idx[0]]

        for it in range(3):
            F1p, F2p = build_F(np.array([[phi1_guess]]), np.array([[phi2_guess]]), coeffs_row, params)
            Fx, Fy = F1p[0, 0], F2p[0, 0]
            Fv = np.array([Fx, Fy])

            h = 0.1
            # derivadas en phi1
            F1p, F2p = build_F(np.array([[phi1_guess + h]]), np.array([[phi2_guess]]), coeffs_row, params) #aca hay error
            F1m, F2m = build_F(np.array([[phi1_guess - h]]), np.array([[phi2_guess]]), coeffs_row, params)
            J11 = (F1p[0, 0] - F1m[0, 0]) / (2 * h)
            J21 = (F2p[0, 0] - F2m[0, 0]) / (2 * h)
            # derivadas en phi2
            F1p, F2p = build_F(np.array([[phi1_guess]]), np.array([[phi2_guess + h]]), coeffs_row, params)
            F1m, F2m = build_F(np.array([[phi1_guess]]), np.array([[phi2_guess - h]]), coeffs_row, params)
            J12 = (F1p[0, 0] - F1m[0, 0]) / (2 * h)
            J22 = (F2p[0, 0] - F2m[0, 0]) / (2 * h)

            J = np.array([[J11, J12], [J21, J22]])

            try:
                delta = np.linalg.solve(J, -Fv)
            except np.linalg.LinAlgError:
                break  # Jacobiano singular -> fold probable

            if np.linalg.norm(delta) > 1.0:
                break  # paso enorme, detener Newton

            phi1_guess += delta[0]
            phi2_guess += delta[1]
            phi1_guess = (phi1_guess + np.pi) % (2 * np.pi) - np.pi
            phi2_guess = (phi2_guess + np.pi) % (2 * np.pi) - np.pi

            if np.linalg.norm(delta) < 1e-5:
                break

        # --- Evaluar estabilidad ---
        try:
            eigvals = np.linalg.eigvals(J)
            re = np.real(eigvals)
            if np.all(re < 0):
                stability = "stable"
            else:
                stability = "saddle"
        except Exception:
            stability = "saddle"

        minima.append((phi1_guess, phi2_guess, d, np.linalg.norm(Fv), stability))
        print(f"Step {s:03d} | d={d:.2f} | |F|={np.linalg.norm(Fv):.2e} | "
              f"phi1={phi1_guess:.3f}, phi2={phi2_guess:.3f} | {stability}")

        d += ds

    minima = np.array(minima, dtype=object)

    # === Gráfico ===
    plt.figure(figsize=(5, 4))
    for phi1_star, phi2_star, d, Fnorm, stab in minima:
        if stab == "stable":
            plt.scatter(d, phi1_star, c='b', s=18, marker='o')
            plt.scatter(d, phi2_star, c='r', s=18, marker='o')
        else:  # saddle
            plt.scatter(d, phi1_star, c='b', s=18, marker='x')
            plt.scatter(d, phi2_star, c='r', s=18, marker='x')

    plt.xlabel(r"$d$")
    plt.ylabel(r"Phase fixed points")
    plt.ylim(-np.pi, np.pi)
    plt.yticks(
        [-np.pi, -np.pi/2, 0, np.pi/2, np.pi],
        [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"]
    )
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()

    return minima



# =========================
#   BLOQUE C — MAIN
# =========================

if __name__ == "__main__":
    directory = r"C:\mnustes_science\simulation_data\FD\PT_dimer\phase_dynamics\alpha=1.000\beta=1.000\mu=0.100\nu=0.200\sigma=3.000\gamma=0.280"  # <--- edita esta línea
    d0 = 21
    mins = continuation(directory, d0=d0, ds=-0.05, steps=400)