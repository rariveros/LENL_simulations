import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
from scipy.special import erf


# Define the real-valued system for integration
def soliton_dynamics_real(t, A_real, omega, V0, sigma, N, d, g):
    A1r, A1i, A2r, A2i = A_real  # Real and imaginary parts

    k = np.sqrt(2 * omega)  # Replace k with sqrt(2*omega)

    # Define erf term for potential interaction
    erf_term = np.sqrt(np.pi) * sigma * erf(k * sigma / 2)

    # Define sinh/cosh interaction term
    sinh_cosh_term = np.sinh(2 * k * d) / np.cosh(2 * k * d) ** 3

    # Compute A1 and A2 in complex form
    A1 = A1r + 1j * A1i
    A2 = A2r + 1j * A2i

    # Define the equations of motion (complex-valued)
    dA1_dt = (k / 2) * (
            - (k ** 2 / 3) * A1 - (4 * N * omega / (3 * N * g)) * A1 ** 3 - (
            4 * N * omega / (3 * N * g)) * A1 * A2 ** 2 * sinh_cosh_term
            + V0 * erf_term * (A1 + A2 * np.exp(-2 * k * d))
    )

    dA2_dt = (k / 2) * (
            - (k ** 2 / 3) * A2 - (4 * N * omega / (3 * N * g)) * A2 ** 3 - (
            4 * N * omega / (3 * N * g)) * A2 * A1 ** 2 * sinh_cosh_term
            + V0 * erf_term * (A2 + A1 * np.exp(-2 * k * d))
    )

    # Convert back to real-valued system
    return [dA1_dt.real, dA1_dt.imag, dA2_dt.real, dA2_dt.imag]


if __name__ == '__main__':

    # Simulation parameters
    omega_values = [1.0]  # Different frequencies
    V0_values = [20]  # Different potential depths
    sigma_values = [1.0]  # Different Gaussian widths
    N_values = [1.0]  # Different interaction strengths
    d_values = [1.0, 5.0, 10.0]  # Different well separations

    t_span = (0, 50)  # Time range
    t_eval = np.linspace(*t_span, 1000)  # Time steps

    # Perform simulations for different parameter sets
    results = []

    for omega in omega_values:
        for V0 in V0_values:
            for sigma in sigma_values:
                for N in N_values:
                    for d in d_values:
                        g = 1.0  # Set interaction parameter g = 1 for simplicity
                        k = np.sqrt(2 * omega)
                        A1_0 = np.sqrt(2 * omega / (N * g))  # Initial condition for A1
                        A2_0 = 0  # Start with A2 = 0

                        # Initial conditions (real and imaginary parts)
                        A_init = [A1_0, 0.0, A2_0, 0.0]

                        # Solve the differential equations
                        sol = solve_ivp(
                            soliton_dynamics_real, t_span, A_init,
                            t_eval=t_eval, args=(omega, V0, sigma, N, d, g), method='RK45'
                        )

                        # Store results
                        results.append((omega, V0, sigma, N, d, sol))

    # Plot results
    plt.figure(figsize=(12, 8))

    for i, (omega, V0, sigma, N, d, sol) in enumerate(results[:6]):  # Plot first few cases
        A1_abs = np.sqrt(sol.y[0] ** 2 + sol.y[1] ** 2)  # Compute |A1|
        A2_abs = np.sqrt(sol.y[2] ** 2 + sol.y[3] ** 2)  # Compute |A2|

        plt.subplot(3, 2, i + 1)
        plt.plot(sol.t, A1_abs, label=r'$|A_1(t)|$')
        plt.plot(sol.t, A2_abs, label=r'$|A_2(t)|$', linestyle='dashed')
        plt.xlabel('Time')
        plt.ylabel('Amplitude')
        plt.title(rf'$\omega={omega}, V_0={V0}, \sigma={sigma}, N={N}, d={d}$')
        plt.legend()

    plt.tight_layout()
    plt.show()