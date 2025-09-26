import matplotlib.pyplot as plt
from functions import *
from back_process import *
from time_integrators import *

from numpy import correlate

if __name__ == '__main__':
    # Hay que barrer g y k, hacer espacio de parametros con fano factor o g.
    # Hay que samplear hartas realizaciones.
    # Hay que tener los diagramas de fase estocasticos y mean-field.
    # Definiendo parámetros

    project_name = '/coherent/test2'
    disc = 'D:/'
    route = 'mnustes_science/simulation_data/FD'
    eq = 'coherent_langevin'
    t_rate = 1

    G2_UU = []
    G2_VV = []
    G2_UV = []

    Fano_UU = []
    Fano_VV = []

    Delta = 0.01
    gamma = 0.1
    Omegas = [0.2]
    K = [0.6] #np.arange(0.5, 0.1, -0.01) # [0.15] #
    g = 0.5

    # Definición de la grilla
    [tmin, tmax, dt] = [0, 4000, 0.1]
    [xmin, xmax, dx] = [0, 10, 1]
    t_grid = np.arange(tmin, tmax + dt, dt)
    x_grid = np.arange(xmin, xmax, dx)
    T = tmax
    Nt = t_grid.shape[0]
    Nx = x_grid.shape[0]
    operators = [0]
    ti = 0.

    sigma = 10.0
    U_init = 0.01 * (np.abs(np.random.normal(0, sigma, Nx)) + 1j * np.abs(np.random.normal(0, sigma, Nx)))
    V_init = 0.01 * (np.abs(np.random.normal(0, sigma, Nx)) + 1j * np.abs(np.random.normal(0, sigma, Nx)))

    U_init = np.concatenate([np.conjugate(U_init), U_init, np.conjugate(U_init), U_init, np.conjugate(U_init), U_init, np.conjugate(U_init), U_init, -np.conjugate(U_init), -U_init, -np.conjugate(U_init), -U_init, -np.conjugate(U_init), -U_init, -np.conjugate(U_init), -U_init])
    V_init = np.concatenate([V_init, np.conjugate(V_init), np.conjugate(V_init), V_init, -V_init, -np.conjugate(V_init), -np.conjugate(V_init), -V_init, V_init, np.conjugate(V_init), np.conjugate(V_init), V_init, -V_init, -np.conjugate(V_init), -np.conjugate(V_init), -V_init,])

    Re1 = np.real(U_init.flatten())
    Im1 = np.imag(U_init.flatten())
    Re2 = np.real(V_init.flatten())
    Im2 = np.imag(V_init.flatten())

    plt.figure(figsize=(6, 5))
    plt.scatter(Re1, Im1, c="k")
    plt.show()
    plt.close()


    for Omega in Omegas:
        for k in K:
            Delta_str = f"{Delta:.{4}f}"
            gamma_str = f"{gamma:.{4}f}"
            Omega_str = f"{Omega:.{4}f}"
            k_str = f"{k:.{4}f}"

            print("####### " + k_str + " #######")

            # Empaquetamiento de parametros, campos y derivadas para integración
            fields_init = [U_init, V_init]
            grids = [t_grid, x_grid, 0]

            parameters_np = np.array([Delta, gamma, Omega, k, g])

            # Midiendo tiempo inicial
            now = datetime.datetime.now()
            print('Hora de Inicio: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
            time_init = time.time()

            final_fields, fields_history, time_grid = RK4_FD(eq, fields_init, parameters_np, grids, dt, Nt, operators, t_rate)

            now = datetime.datetime.now()
            print('Hora de Término: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
            time_fin = time.time()
            print(str(time_fin - time_init) + ' seg')

            # Reobteniendo campos
            U = np.array(fields_history)[:, 0]
            V = np.array(fields_history)[:, 1]

            # Guardando datos
            file = disc + route + project_name
            subfile = "/Delta=" + Delta_str + "/gamma=" + gamma_str + "/Omega=" + Omega_str + "/k=" + k_str
            if not os.path.exists(file + subfile):
                os.makedirs(file + subfile)

            np.savetxt(file + subfile + '/U.txt', U, delimiter=',')
            np.savetxt(file + subfile + '/V.txt', V, delimiter=',')
            np.savetxt(file + subfile + '/parameters.txt', parameters_np, delimiter=',')
            np.savetxt(file + subfile + '/T.txt', t_grid, delimiter=',')

            lightness = 1
            U_light = U[0::lightness]
            V_light = V[0::lightness]
            t_light = time_grid[0::lightness]

            # --- eliminar 20% inicial (transiente) ---
            cut = int(0.2 * len(t_light))
            U_light = U_light[cut:]
            V_light = V_light[cut:]
            t_light = t_light[cut:]

            # FFT along time axis
            spec_U = np.fft.fftshift(np.fft.fft(U_light, axis=0), axes=0)
            spec_V = np.fft.fftshift(np.fft.fft(V_light, axis=0), axes=0)

            # Power spectrum per realization
            psd_U = np.abs(spec_U) ** 2
            psd_V = np.abs(spec_V) ** 2

            # Average over realizations (columns)
            avg_psd_U = np.mean(psd_U, axis=1)
            avg_psd_V = np.mean(psd_V, axis=1)

            # Frequency grid
            freqs = np.fft.fftshift(np.fft.fftfreq(len(t_light), d=dt))

            plt.figure()
            plt.plot(freqs, avg_psd_U, 'b', label='U avg')
            plt.plot(freqs, avg_psd_V, 'r', label='V avg')
            plt.xlabel("Frequency")
            plt.ylabel("Average Power Spectrum")
            plt.grid(linestyle='--', alpha=0.5)
            plt.xlim(-0.5, 0.5)
            plt.legend()
            plt.savefig("spectrum_avg.png", dpi=200)
            plt.close()

            """
            # === Time series (con transiente eliminado) ===
            fig, (ax1, ax2) = plt.subplots(2)
            ax1.plot(t_light, np.real(U_light), c="b", label="$\\textrm{Re}\\, \\alpha_i$")
            ax1.plot(t_light, np.imag(U_light), c="r", label="$\\textrm{Im}\\, \\alpha_i$")
            ax1.set_xlabel('$t$', size=20)
            ax1.set_ylabel('$\\alpha_1$', size=20)
            ax1.set_xlim([t_light[0], t_light[-1]])
            ax1.grid(linestyle='--', alpha=0.5)
            ax1.legend(loc="upper right", fontsize=12)

            ax2.plot(t_light, np.real(V_light), c="b")
            ax2.plot(t_light, np.imag(V_light), c="r")
            ax2.set_xlabel('$t$', size=20)
            ax2.set_ylabel('$\\alpha_2$', size=20)
            ax2.set_xlim([t_light[0], t_light[-1]])
            ax2.grid(linestyle='--', alpha=0.5)
            plt.savefig(file + subfile + "/timeseries.png", dpi=300)
            plt.close()

            # === Photon number statistics (post-transient) ===
            nU = np.abs(U_light.flatten()) ** 2
            nV = np.abs(V_light.flatten()) ** 2
            mean_nU = np.mean(nU)
            mean_nV = np.mean(nV)

            # === Fano factors ===
            fano_U = (np.mean(nU ** 2) - mean_nU ** 2) / mean_nU
            fano_V = (np.mean(nV ** 2) - mean_nV ** 2) / mean_nV

            Fano_UU.append(fano_U)
            Fano_VV.append(fano_V)

            N = len(nU)

            taus = (t_light - t_light[0])  # tau directamente desde el grid de tiempo

            # === g2_UU ===
            acorr_U = np.correlate(nU, nU, mode='full')[N - 1:]
            g2_UU = acorr_U / ((N - np.arange(N)) * mean_nU ** 2)

            plt.figure()
            plt.plot(taus[:len(g2_UU)], g2_UU, 'k')
            plt.xlabel('$\\tau$')
            plt.ylabel('$g^{2}_{UU}(\\tau)$')
            plt.grid(linestyle='--', alpha=0.5)
            plt.xlim(0, 400)
            plt.savefig(file + subfile + "/g2_UU.png", dpi=150)
            plt.close()

            # === g2_VV ===
            acorr_V = np.correlate(nV, nV, mode='full')[N - 1:]
            g2_VV = acorr_V / ((N - np.arange(N)) * mean_nV ** 2)

            plt.figure()
            plt.plot(taus[:len(g2_VV)], g2_VV, 'k')
            plt.xlabel('$\\tau$')
            plt.ylabel('$g^{2}_{VV}(\\tau)$')
            plt.grid(linestyle='--', alpha=0.5)
            plt.xlim(0, 400)
            plt.savefig(file + subfile + "/g2_VV.png", dpi=150)
            plt.close()

            # === g2_UV (cross correlation) ===
            cross_UV = np.correlate(nU, nV, mode='full')[N - 1:]
            g2_UV = cross_UV / ((N - np.arange(N)) * mean_nU * mean_nV)

            plt.figure()
            plt.plot(taus[:len(g2_UV)], g2_UV, 'k')
            plt.xlabel('$\\tau$')
            plt.ylabel('$g^{2}_{UV}(\\tau)$')
            plt.grid(linestyle='--', alpha=0.5)
            plt.xlim(0, 400)
            plt.savefig(file + subfile + "/g2_UV.png", dpi=150)
            plt.close()
            """
            Re1 = np.real(U_light.flatten())
            Im1 = np.imag(U_light.flatten())
            Re2 = np.real(V_light.flatten())
            Im2 = np.imag(V_light.flatten())

            print(Re1)

            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4, 8))
            hist_01 = ax1.hist2d(Re1, Im2, bins=160, density=True, cmap="inferno")
            ax1.set_xlabel("$\\textrm{Re }(\\alpha_1)$")
            ax1.set_ylabel("$\\textrm{Im }(\\alpha_1)$")

            hist_01 = ax2.hist2d(Re2, Im1, bins=160, density=True, cmap="inferno")
            ax2.set_xlabel("$\\textrm{Re }(\\alpha_2)$")
            ax2.set_ylabel("$\\textrm{Im }(\\alpha_2)$")
            #plt.gca().set_aspect("equal")
            plt.savefig("phase_space.png", dpi=300)
            plt.close()

            #g2_0_UU = np.mean(nU ** 2) / (mean_nU ** 2)
            #g2_0_VV = np.mean(nV ** 2) / (mean_nV ** 2)
            #g2_0_UV = np.mean(nU * nV) / (mean_nU * mean_nV)

            #G2_UU.append(g2_0_UU)
            #G2_VV.append(g2_0_VV)
            #G2_UV.append(g2_0_UV)

            U_init = U_light[-1] + 0.01 * (np.random.rand(16 * Nx) + 1j * np.random.rand(16 * Nx))
            V_init = V_light[-1] + 0.01 * (np.random.rand(16 * Nx) + 1j * np.random.rand(16 * Nx))
    """
    plt.figure()
    plt.scatter(K, G2_UU, c="b", label="$UU$")
    plt.scatter(K, G2_VV, c="r", label="$UU$")
    plt.scatter(K, G2_UV, c="k", label="$UU$")
    plt.xlabel('$\kappa$')
    plt.ylabel('$g^{2}(0)$')
    plt.grid(linestyle='--', alpha=0.5)
    plt.savefig(file + "/g2.png", dpi=150)
    plt.close()

    plt.figure()
    plt.plot(K, Fano_UU, 'bo-', label='U')
    plt.plot(K, Fano_VV, 'ro-', label='V')
    plt.xlabel('$\kappa$')
    plt.ylabel('$\\textrm{Fano factor}$')
    plt.grid(linestyle='--', alpha=0.5)
    plt.legend()
    plt.savefig(file + "/fano.png", dpi=150)
    plt.close()
    """
