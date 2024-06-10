from back_process import *


if __name__ == '__main__':
    directory = 'C:/mnustes_science/simulation_data/FD/PT_extended/variando_distancia_rabi/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.0180/sigma=6.000/gamma=0.1850/dist=30.2920'

    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    X = np.loadtxt(directory + '/X.txt', delimiter=',')
    T = np.loadtxt(directory + '/T.txt', delimiter=',')
    params = np.loadtxt(directory + '/parameters.txt', delimiter=',')
    # [alpha, beta, gamma_0, dist,  mu, nu, sigma]
    dist = params[3]
    beta = 0.004811649356064012

    Z_r = Z_r / np.sqrt(beta)
    Z_i = Z_i / np.sqrt(beta)

    Z = Z_r + 1j * Z_i
    Z_conj = Z_r - 1j * Z_i
    Z_modulo = np.absolute(Z)
    center_L = - 0.5 * params[3]
    center_R = 0.5 * params[3]
    L = X[-1] - X[0]

    Nx = len(X)
    Nt = len(T)

    ti, tf = int(0.1 * Nt), Nt
    Z_r = Z_r[ti:tf, :]
    Z_i = Z_i[ti:tf, :]
    Z_modulo = Z_modulo[ti:tf, :]
    T = T[ti:tf] - T[ti]
    dt = T[1] - T[0]

    Nx = len(X)
    Nt = len(T)

    ###      CORRELATION AND PERIOD     ###

    Z_mod_L = Z_modulo[:, 0:int(Nx / 2)]
    Z_mod_R = Z_modulo[:, int(Nx / 2) + 1:]
    X_L = X[0:int(Nx / 2)]
    X_R = X[int(Nx / 2) + 1:]

    N_L = integrate.simpson(Z_mod_L ** 2, X_L) / L
    N_R = integrate.simpson(Z_mod_R ** 2, X_R) / L

    CCF = np.correlate(N_L - np.mean(N_L), N_R - np.mean(N_R), "full")
    tau = np.arange(-T[-1], T[-1], dt)
    Ntau = len(tau)
    dtau = tau[1] - tau[0]

    CCF_max, tau_max, CCF_I = max_finder(CCF[:-1], tau, Ntau, dtau)
    mod_strobo = []
    real_strobo = []
    imag_strobo = []
    tau_R = []
    maxval = np.amax(CCF)
    for i in range(len(tau_max)):
        if CCF_max[i] > 0.25 * maxval:
            tau_R.append(tau_max[i])
    period = np.mean(np.diff(tau_R))
    period_std = np.std(np.diff(tau_R))

    fig, ((ax1), (ax2)) = plt.subplots(nrows=2, ncols=1)
    I_period = np.argmin(np.abs(T - period))
    I_initial = np.argmax(N_L[:int(2 * I_period)])
    t_initial = T[I_initial]
    t_actual = [t_initial]
    n = 0

    while t_actual[-1] + 2 * period < T[-1]:
        t_actual.append(t_actual[-1] + 1 * period)
        i_actual = np.argmin(np.abs(T - (t_actual[-1] + 1 * period)))
        mod_strobo.append(Z_modulo[i_actual, :])
        real_strobo.append(Z_r[i_actual, :])
        imag_strobo.append(Z_i[i_actual, :])

        if n % 2 == 0:
            ls = "--"
        else:
            ls = "solid"
        if n == 0:
            #ax1.plot(X, Z_r[i_actual, :], color="b", linestyle=ls, label="$\psi_I$")
            #ax1.plot(X, Z_i[i_actual, :], color="r", linestyle=ls, label="$\psi_R$")
            ax1.plot(X, Z_modulo[i_actual, :], color="k", label="$|\psi|$")

        else:
            #ax1.plot(X, Z_r[i_actual, :], color="b", linestyle=ls)
            #ax1.plot(X, Z_i[i_actual, :], color="r", linestyle=ls)
            ax1.plot(X, Z_modulo[i_actual, :], color="k")
        n = n + 1

    I_initial = np.argmax(N_R[:int(2 * I_period)])
    t_initial = T[I_initial]
    t_actual = [t_initial]
    n = 0
    while t_actual[-1] + 2 * period < T[-1]:
        t_actual.append(t_actual[-1] + 1 * period)
        i_actual = np.argmin(np.abs(T - (t_actual[-1] + 1 * period)))
        mod_strobo.append(Z_modulo[i_actual, :])
        real_strobo.append(Z_r[i_actual, :])
        imag_strobo.append(Z_i[i_actual, :])

        if n % 2 == 0:
            ls = "--"
        else:
            ls = "solid"
        #ax2.plot(X, Z_r[i_actual, :], color="b", linestyle=ls)
        #ax2.plot(X, Z_i[i_actual, :], color="r", linestyle=ls)
        ax1.plot(X, Z_modulo[i_actual, :], color="k")
        n = n + 1
    print(len(mod_strobo))
    ax1.text(-140, 3, "$d=" + str(dist) + "$", fontsize=15, rotation=0, zorder=13)
    ax1.tick_params(axis="x", direction="in", labeltop=False, labelbottom=False, top=True, bottom=True, grid_markersize=15)
    ax1.set_xlim(X[0], X[-1])
    ax1.grid(alpha=0.3)
    ax1.legend(fontsize=13, loc="upper right")
    np.savetxt("mod_strobo_rabi.txt", mod_strobo, delimiter=',')
    ax2.tick_params(axis="x", direction="in", labeltop=False, labelbottom=True, top=True, bottom=True, grid_markersize=15)
    ax2.set_xlim(X[0], X[-1])
    ax2.set_xlabel("$x$", fontsize=25)
    ax2.grid(alpha=0.3)
    plt.savefig("strobo_" + str(dist) + ".png", dpi=250)