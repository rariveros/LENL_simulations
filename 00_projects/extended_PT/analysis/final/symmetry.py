import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    principal_sim = "C:/mnustes_science/simulation_data/FD/PDNLS_extended_PT/extras"
    dir_02 = principal_sim + "/dimensional_dist_02/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.0180/sigma=6.000/gamma=0.1850/dist=30.2930"
    dir_01 = principal_sim + "/dimensional_dist_02/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.0180/sigma=6.000/gamma=0.1850/dist=22.0000"

    Z_r = np.loadtxt(dir_01 + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(dir_01 + '/field_img.txt', delimiter=',')
    X = np.loadtxt(dir_01 + '/X.txt', delimiter=',')
    T = np.loadtxt(dir_01 + '/T.txt', delimiter=',')

    Nt = len(T)
    Nx = len(X)
    dt = T[1] - T[0]

    Z_r = Z_r[int(Nt/2):, :]
    Z_i = Z_i[int(Nt / 2):, :]
    T = T[int(Nt / 2):]
    T = T - T[0]
    Z_complex = Z_r + 1j * Z_i
    Z_modulo = np.abs(Z_complex)

    ###

    N_R = integrate.simpson(np.abs(Z_r) ** 2, X)

    CCF = np.correlate(N_R, N_R, "full")

    tau = np.arange(-T[-1], T[-1] + dt, dt)
    Ntau = len(tau)
    dtau = tau[1] - tau[0]

    CCF_max, tau_max, CCF_I = max_finder(CCF, tau, Ntau, dtau)

    tau_R = []
    maxval = np.amax(CCF)
    for i in range(len(tau_max)):
        if CCF_max[i] > 0.8 * maxval:
            tau_R.append(tau_max[i])
    period = np.mean(np.diff(tau_R))
    period_std = np.std(np.diff(tau_R))
    t_dist = np.abs(period - T)
    I_period = np.argmin(t_dist)
    print(I_period)

    Nt = len(T)
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)

    pc_01 = ax1.pcolor(X, T, Z_r, cmap="viridis")
    pc_02 = ax2.pcolor(X, T-period/2, -np.flip(Z_i, axis=1), cmap="plasma")
    pc_03 = ax3.pcolor(X, T, Z_i, cmap="viridis")
    pc_04 = ax4.pcolor(X, T-period/2, np.flip(Z_r, axis=1), cmap="plasma")
    ax1.set_ylim([0, 400])
    ax2.set_ylim([0, 400])
    ax3.set_ylim([0, 400])
    ax4.set_ylim([0, 400])
    plt.show()
    plt.close()
