import numpy as np
import matplotlib.pyplot as plt

from scipy.special import factorial

if __name__ == '__main__':
    # PARAMETERS
    N = 10
    Delta = 0.01
    Omega = 0.01
    kappa = 0.2
    g = 0.01
    gamma = 0.001
    T = 3000
    dt = 0.1
    steps = int(T / dt)
    times = np.linspace(0, T, steps)
    plot_start = 0 #-int(500 / dt)


    # OPERATORS
    def create_annihilation_op(N):
        a = np.zeros((N, N), dtype=complex)
        for n in range(1, N):
            a[n - 1, n] = np.sqrt(n)
        return a


    a = create_annihilation_op(N)
    adag = a.conj().T
    I = np.eye(N)
    b1 = np.kron(a, I)
    b2 = np.kron(I, a)
    b1_dag = b1.conj().T
    b2_dag = b2.conj().T
    b1_sq = b1 @ b1
    b2_sq = b2 @ b2
    b1_dag_sq = b1_dag @ b1_dag
    b2_dag_sq = b2_dag @ b2_dag
    n1_op = b1_dag @ b1
    n2_op = b2_dag @ b2


    # INITIAL STATE: coherent with π/2 phase difference
    def coherent_state(alpha, N):
        state = sum(
            (alpha ** n / np.sqrt(factorial(n))) * np.eye(N)[:, n:n + 1]
            for n in range(N)
        )
        return state / np.linalg.norm(state)

    phase = 0.5 * np.pi
    alpha1 = 0.0 + 0.0 * 1j
    alpha2 = 0.1 * np.exp(-1j * phase) # alpha1 * 1j #
    coh1 = coherent_state(alpha1, N)
    coh2 = coherent_state(alpha2, N)
    psi0 = np.kron(coh1, coh2)
    rho0 = psi0 @ psi0.conj().T

    #vac = np.zeros((N, 1), dtype=complex)
    #vac[0, 0] = 1
    #one_photon = np.zeros((N, 1), dtype=complex)
    #one_photon[1, 0] = 1

    #psi0 = np.kron(one_photon, vac)
    #rho0 = psi0 @ psi0.conj().T

    # HAMILTONIAN
    H = (
            Delta * (b1_dag @ b1 + b2_dag @ b2)
            + Omega * (b1_sq + b1_dag_sq - b2_sq - b2_dag_sq)
            - kappa * (b1 @ b2_dag + b1_dag @ b2)
            + g * (b1_dag_sq @ b1_sq + b2_dag_sq @ b2_sq)
    )

    # LINDLBLAD TERMS
    L1 = np.sqrt(gamma) * b1
    L2 = np.sqrt(gamma) * b2


    def lindblad_rhs(rho):
        comm = -1j * (H @ rho - rho @ H)
        D1 = L1 @ rho @ L1.conj().T - 0.5 * (L1.conj().T @ L1 @ rho + rho @ L1.conj().T @ L1)
        D2 = L2 @ rho @ L2.conj().T - 0.5 * (L2.conj().T @ L2 @ rho + rho @ L2.conj().T @ L2)
        return comm + D1 + D2


    def rk4_step(rho, dt):
        k1 = lindblad_rhs(rho)
        k2 = lindblad_rhs(rho + 0.5 * dt * k1)
        k3 = lindblad_rhs(rho + 0.5 * dt * k2)
        k4 = lindblad_rhs(rho + dt * k3)
        return rho + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)


    # SIMULATION
    rho = rho0.copy()
    n1_vals = np.zeros(steps)
    n2_vals = np.zeros(steps)
    b1_vals = np.zeros(steps, dtype=complex)
    b2_vals = np.zeros(steps, dtype=complex)

    for i in range(steps):
        n1_vals[i] = np.real(np.trace(rho @ n1_op))
        n2_vals[i] = np.real(np.trace(rho @ n2_op))
        b1_vals[i] = np.trace(rho @ b1)
        b2_vals[i] = np.trace(rho @ b2)
        rho = rk4_step(rho, dt)

    # EXTRACT LAST 500 UNITS

    times_tail = times[plot_start:]
    n1_tail = n1_vals[plot_start:]
    n2_tail = n2_vals[plot_start:]
    b1_tail = b1_vals[plot_start:]
    b2_tail = b2_vals[plot_start:]
    re_b1_tail = np.real(b1_tail)
    im_b1_tail = np.imag(b1_tail)
    re_b2_tail = np.real(b2_tail)
    im_b2_tail = np.imag(b2_tail)
    phase_diff_tail = (np.angle(b1_tail) - np.angle(b2_tail) + np.pi) % (2 * np.pi) - np.pi

    # PLOT
    plt.figure(figsize=(14, 8))

    plt.subplot(2, 2, 1)
    plt.plot(times_tail, n1_tail, label="⟨n₁⟩", color="tab:orange")
    plt.plot(times_tail, n2_tail, label="⟨n₂⟩", color="tab:red")
    plt.title("Populations (final 500 units)")
    plt.xlabel("Time")
    plt.ylabel("⟨n⟩")
    plt.legend()
    plt.grid(True)

    plt.subplot(2, 2, 2)
    plt.plot(times_tail, re_b1_tail, label="Re⟨b₁⟩", color="tab:blue")
    plt.plot(times_tail, re_b2_tail, label="Re⟨b₂⟩", color="tab:green")
    plt.title("Re⟨b⟩ (final 500 units)")
    plt.xlabel("Time")
    plt.ylabel("Re⟨b⟩")
    plt.legend()
    plt.grid(True)

    plt.subplot(2, 2, 3)
    plt.plot(times_tail, im_b1_tail, label="Im⟨b₁⟩", color="tab:blue")
    plt.plot(times_tail, im_b2_tail, label="Im⟨b₂⟩", color="tab:green")
    plt.title("Im⟨b⟩ (final 500 units)")
    plt.xlabel("Time")
    plt.ylabel("Im⟨b⟩")
    plt.legend()
    plt.grid(True)

    plt.subplot(2, 2, 4)
    plt.plot(times_tail, phase_diff_tail, label="arg⟨b₁⟩ - arg⟨b₂⟩", color='tab:purple')
    plt.axhline(np.pi / 2, color='gray', linestyle='--', linewidth=1, label="π/2")
    plt.axhline(0, color='black', linestyle=':', linewidth=1)
    plt.title("Phase Difference (final 500 units)")
    plt.xlabel("Time")
    plt.ylabel("Δ Phase (rad)")
    plt.legend()
    plt.grid(True)

    plt.suptitle("Final 500 Time Units (T = 3000, Initial π/2 Phase Difference)")
    plt.tight_layout()
    plt.show()