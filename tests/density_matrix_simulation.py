import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    # Pauli matrices and lowering operator
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
    sigma_p = np.array([[0, 1], [0, 0]], dtype=complex)
    sigma_m = np.array([[0, 0], [1, 0]], dtype=complex)

    # Parameters
    omega = 1.0     # Hamiltonian frequency
    gamma = 0.2     # Decay rate

    # Hamiltonian and Lindblad operator
    H = 0.5 * omega * sigma_z
    L = np.sqrt(gamma) * sigma_m

    # Initial state: |+x⟩ = (|0⟩ + |1⟩)/√2
    psi0 = (1/np.sqrt(2)) * np.array([[1], [1]], dtype=complex)
    rho0 = psi0 @ psi0.conj().T

    # Lindblad RHS
    def lindblad_rhs(rho):
        comm = -1j * (H @ rho - rho @ H)
        dissip = L @ rho @ L.conj().T - 0.5 * (L.conj().T @ L @ rho + rho @ L.conj().T @ L)
        return comm + dissip

    # RK4 step for matrix-valued ODE
    def rk4_step(rho, dt):
        k1 = lindblad_rhs(rho)
        k2 = lindblad_rhs(rho + 0.5 * dt * k1)
        k3 = lindblad_rhs(rho + 0.5 * dt * k2)
        k4 = lindblad_rhs(rho + dt * k3)
        return rho + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4)

    # Time evolution
    T = 20
    dt = 0.01
    steps = int(T / dt)
    times = np.linspace(0, T, steps)

    # Arrays to store expectation values
    expect_x = np.zeros(steps)
    expect_y = np.zeros(steps)
    expect_z = np.zeros(steps)

    rho = rho0.copy()
    for i in range(steps):
        expect_x[i] = np.real(np.trace(rho @ sigma_x))
        expect_y[i] = np.real(np.trace(rho @ sigma_y))
        expect_z[i] = np.real(np.trace(rho @ sigma_z))
        rho = rk4_step(rho, dt)

    # Plot
    plt.figure(figsize=(8, 4))
    plt.plot(times, expect_x, label='⟨σx⟩')
    plt.plot(times, expect_y, label='⟨σy⟩')
    plt.plot(times, expect_z, label='⟨σz⟩')
    plt.xlabel("Time")
    plt.ylabel("Expectation values")
    plt.legend()
    plt.grid(True)
    plt.title("Spin-1/2 Lindblad Evolution (Decay)")
    plt.tight_layout()
    plt.show()