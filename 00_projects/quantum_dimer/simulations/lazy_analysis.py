import matplotlib.pyplot as plt
import numpy as np
import os, time, datetime

from functions import *
from back_process import *
from time_integrators import *
from numpy import correlate

if __name__ == '__main__':
    # Parameters
    project_name = '/coherent/test2'
    disc = 'D:/'
    route = 'mnustes_science/simulation_data/FD'
    eq = 'coherent_langevin'
    t_rate = 1

    Delta = 0.01
    gamma = 0.1
    Omegas = [0.2]
    K = np.arange(0, 0.61, 0.1) #[0.6]  # or np.arange(...) if you want to sweep
    g = 0.1

    # Grids
    [tmin, tmax, dt] = [0, 2000, 0.1]
    [xmin, xmax, dx] = [0, 4, 1]
    t_grid = np.arange(tmin, tmax + dt, dt)
    x_grid = np.arange(xmin, xmax, dx)
    T = tmax
    Nt = t_grid.shape[0]
    Nx = x_grid.shape[0]
    operators = [0]
    ti = 0.

    sigma = 10.0
    U_init = 0.01 * (np.abs(np.random.normal(0, sigma, Nx)) +
                     1j * np.abs(np.random.normal(0, sigma, Nx)))
    V_init = 0.01 * (np.abs(np.random.normal(0, sigma, Nx)) +
                     1j * np.abs(np.random.normal(0, sigma, Nx)))

    U_init = np.concatenate([np.conjugate(U_init), U_init,
                             np.conjugate(U_init), U_init,
                             np.conjugate(U_init), U_init,
                             np.conjugate(U_init), U_init,
                             -np.conjugate(U_init), -U_init,
                             -np.conjugate(U_init), -U_init,
                             -np.conjugate(U_init), -U_init,
                             -np.conjugate(U_init), -U_init])
    V_init = np.concatenate([V_init, np.conjugate(V_init),
                             np.conjugate(V_init), V_init,
                             -V_init, -np.conjugate(V_init),
                             -np.conjugate(V_init), -V_init,
                             V_init, np.conjugate(V_init),
                             np.conjugate(V_init), V_init,
                             -V_init, -np.conjugate(V_init),
                             -np.conjugate(V_init), -V_init])

    # === lists for final scatter and cmap ===
    AvgOcc_U, AvgOcc_V = [], []
    VarOcc_U, VarOcc_V = [], []
    Ks = []

    Spectra_U, Spectra_V = [], []
    Freqs = None

    for Omega in Omegas:
        for k in K:
            Delta_str = f"{Delta:.{4}f}"
            gamma_str = f"{gamma:.{4}f}"
            Omega_str = f"{Omega:.{4}f}"
            k_str = f"{k:.{4}f}"

            print("####### " + k_str + " #######")

            fields_init = [U_init, V_init]
            grids = [t_grid, x_grid, 0]
            parameters_np = np.array([Delta, gamma, Omega, k, g])

            now = datetime.datetime.now()
            print('Hora de Inicio: ' + str(now.hour) + ':' +
                  str(now.minute) + ':' + str(now.second))
            time_init = time.time()

            final_fields, fields_history, time_grid = RK4_FD(
                eq, fields_init, parameters_np, grids, dt, Nt, operators, t_rate
            )

            now = datetime.datetime.now()
            print('Hora de Término: ' + str(now.hour) + ':' +
                  str(now.minute) + ':' + str(now.second))
            time_fin = time.time()
            print(str(time_fin - time_init) + ' seg')

            # Reobteniendo campos
            U = np.array(fields_history)[:, 0]
            V = np.array(fields_history)[:, 1]

            # Guardar datos
            file = disc + route + project_name
            subfile = "/Delta=" + Delta_str + "/gamma=" + gamma_str + "/Omega=" + Omega_str + "/k=" + k_str
            if not os.path.exists(file + subfile):
                os.makedirs(file + subfile)

            np.savetxt(file + subfile + '/U.txt', U, delimiter=',')
            np.savetxt(file + subfile + '/V.txt', V, delimiter=',')
            np.savetxt(file + subfile + '/parameters.txt',
                       parameters_np, delimiter=',')
            np.savetxt(file + subfile + '/T.txt', t_grid, delimiter=',')

            # === cut transient ===
            U_light = U
            V_light = V
            t_light = time_grid
            cut = int(0.2 * len(t_light))
            U_light = U_light[cut:]
            V_light = V_light[cut:]
            t_light = t_light[cut:]

            # === order parameters per realization ===
            for j in range(U_light.shape[1]):  # each realization
                nU = np.abs(U_light[:, j])**2
                nV = np.abs(V_light[:, j])**2

                # averages
                AvgOcc_U.append(np.mean(nU))
                AvgOcc_V.append(np.mean(nV))

                # variances
                VarOcc_U.append(np.var(nU))
                VarOcc_V.append(np.var(nV))

                # store k
                Ks.append(k)

                # spectra
                spec_U = np.fft.fftshift(np.fft.fft(U_light[:, j]))
                spec_V = np.fft.fftshift(np.fft.fft(V_light[:, j]))
                psd_U = np.abs(spec_U)**2
                psd_V = np.abs(spec_V)**2

                # normalize
                psd_U /= np.max(psd_U)
                psd_V /= np.max(psd_V)

                Spectra_U.append(psd_U)
                Spectra_V.append(psd_V)

                if Freqs is None:
                    Freqs = np.fft.fftshift(
                        np.fft.fftfreq(len(t_light), d=dt))

    # === Example plots ===

    # Scatter plot: average occupation vs coupling
    plt.figure()
    plt.scatter(Ks, AvgOcc_U, c='b', label='U avg occ')
    plt.scatter(Ks, AvgOcc_V, c='r', label='V avg occ')
    plt.xlabel('$\\kappa$')
    plt.ylabel('Average occupation')
    plt.legend()
    plt.grid(linestyle='--', alpha=0.5)
    plt.savefig("avg_occ_vs_k.png", dpi=200)
    plt.close()

    # Scatter plot: variance vs coupling
    plt.figure()
    plt.scatter(Ks, VarOcc_U, c='b', label='U var')
    plt.scatter(Ks, VarOcc_V, c='r', label='V var')
    plt.xlabel('$\\kappa$')
    plt.ylabel('Variance of intensity')
    plt.legend()
    plt.grid(linestyle='--', alpha=0.5)
    plt.savefig("var_vs_k.png", dpi=200)
    plt.close()

    # Cmap of spectrum vs kappa (example with U)
    Spectra_U = np.array(Spectra_U)  # shape (Nreal_total, Nfreq)
    Spectra_V = np.array(Spectra_V)
    Ks_arr = np.array(Ks)

    plt.figure(figsize=(6, 5))
    plt.pcolormesh(Freqs, Ks_arr, Spectra_U,
                   shading='auto', cmap='inferno')
    plt.xlabel("Frequency")
    plt.ylabel("$\\kappa$")
    plt.colorbar(label="$S_U(\omega)$")
    plt.savefig("spectrum_vs_k_U.png", dpi=200)
    plt.close()

    plt.figure(figsize=(6, 5))
    plt.pcolormesh(Freqs, Ks_arr, Spectra_V,
                   shading='auto', cmap='inferno')
    plt.xlabel("Frequency")
    plt.ylabel("$\\kappa$")
    plt.colorbar(label="$S_V(\omega)$")
    plt.savefig("spectrum_vs_k_V.png", dpi=200)
    plt.close()
