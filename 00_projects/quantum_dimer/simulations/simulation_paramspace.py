import matplotlib.pyplot as plt
import numpy as np
import os, time, datetime
from functions import *
from back_process import *
from time_integrators import *
from numpy import correlate
from matplotlib.colors import LogNorm
from scipy.ndimage import gaussian_filter1d  # <-- para suavizar espectros

if __name__ == '__main__':
    # Parameters
    project_name = '/coherent/parameter_space_01'
    disc = 'C:/'
    route = 'mnustes_science/simulation_data/FD'
    eq = 'coherent_langevin'
    t_rate = 1

    Delta = 0.0
    gamma = 0.1
    Omega = 0.2
    K = np.arange(0.0, 1.0, 0.01)
    #gs = [0.0011]

    # Grids
    [tmin, tmax, dt] = [0, 10000, 0.1]
    [xmin, xmax, dx] = [0, 10, 1]
    t_grid = np.arange(tmin, tmax + dt, dt)
    x_grid = np.arange(xmin, xmax, dx)
    Nt = t_grid.shape[0]
    Nx = x_grid.shape[0]
    operators = [0]

    sigma = 10.0
    U_init = 0.01 * (np.abs(np.random.normal(0, sigma, Nx)) + 1j * np.abs(np.random.normal(0, sigma, Nx)))
    V_init = 0.01 * (np.abs(np.random.normal(0, sigma, Nx)) + 1j * np.abs(np.random.normal(0, sigma, Nx)))

    U_init = np.concatenate([np.conjugate(U_init), U_init, np.conjugate(U_init), U_init, np.conjugate(U_init), U_init, np.conjugate(U_init), U_init, -np.conjugate(U_init), -U_init, -np.conjugate(U_init), -U_init, -np.conjugate(U_init), -U_init, -np.conjugate(U_init), -U_init])
    V_init = np.concatenate([V_init, np.conjugate(V_init), np.conjugate(V_init), V_init, -V_init, -np.conjugate(V_init), -np.conjugate(V_init), -V_init, V_init, np.conjugate(V_init), np.conjugate(V_init), V_init, -V_init, -np.conjugate(V_init), -np.conjugate(V_init), -V_init,])

    QUANTUMNESS = np.arange(2.0, -1.05, -0.05)# [2.0, 1.0, 0.0, -1.0] #
    output_data = []
    for quantumness in QUANTUMNESS:
        g = gamma * 10 ** (-quantumness)
        AvgOcc_U, AvgOcc_V = [], []
        VarOcc_U, VarOcc_V = [], []
        Ks = []
        Ks_mean = []

        Spectra_U, Spectra_V = [], []
        for k in K:
            Delta_str = f"{Delta:.{4}f}"
            gamma_str = f"{gamma:.{4}f}"
            Omega_str = f"{Omega:.{4}f}"
            k_str = f"{k:.{4}f}"
            g_str = f"{g:.{4}f}"

            print("####### κ = ", k, " #######")
            # Midiendo tiempo inicial
            now = datetime.datetime.now()
            print('Hora de Inicio: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
            time_init = time.time()

            fields_init = [U_init, V_init]
            grids = [t_grid, x_grid, 0]
            parameters_np = np.array([Delta, gamma, Omega, k, g])

            final_fields, fields_history, time_grid = RK4_FD(eq, fields_init, parameters_np, grids, dt, Nt, operators, t_rate)

            # Reobteniendo campos
            U = np.array(fields_history)[:, 0]
            V = np.array(fields_history)[:, 1]

            # cut transient
            cut = int(0.5 * len(time_grid))
            U_light = U[cut:]
            V_light = V[cut:]
            t_light = time_grid[cut:]

            # === acumular espectros brutos por cada κ ===
            psd_U_all = []
            psd_V_all = []

            for j in range(U_light.shape[1]):
                nU = np.abs(U_light[:, j]) ** 2
                nV = np.abs(V_light[:, j]) ** 2

                # order parameters
                AvgOcc_U.append(np.mean(nU))
                AvgOcc_V.append(np.mean(nV))
                VarOcc_U.append(np.var(nU))
                VarOcc_V.append(np.var(nV))
                Ks.append(k)

                # FFT con zero-padding
                Npad = 4 * len(t_light)
                spec_U = np.fft.fftshift(np.fft.fft(U_light[:, j], n=Npad))
                spec_V = np.fft.fftshift(np.fft.fft(V_light[:, j], n=Npad))

                psd_U = np.abs(spec_U) ** 2
                psd_V = np.abs(spec_V) ** 2

                psd_U_all.append(psd_U)
                psd_V_all.append(psd_V)

            # === promedio sobre realizaciones ===
            psd_U_mean = np.mean(psd_U_all, axis=0)
            psd_V_mean = np.mean(psd_V_all, axis=0)

            Freqs = np.fft.fftshift(np.fft.fftfreq(Npad, d=dt))

            spec = psd_U_mean

            idx_peak = np.argmax(spec)
            f_peak = Freqs[idx_peak]
            S_max = spec[idx_peak]
            S_noise = spec.mean()
            M = S_max

            # === normalizar después del promedio (cada κ con máximo = 1) ===
            psd_U_mean /= np.max(psd_U_mean)

            # === suavizado ligero ===
            #psd_U_mean = gaussian_filter1d(psd_U_mean, sigma=10)
            #psd_V_mean = gaussian_filter1d(psd_V_mean, sigma=5)
            del U_light, V_light, U, V
            output_data.append([k, quantumness, g, f_peak, M, S_noise])

            now = datetime.datetime.now()
            print('Hora de Término: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
            time_fin = time.time()
            print(str(time_fin - time_init) + ' seg')
        file = disc + route + project_name
        subfile = "/Delta=" + Delta_str + "/gamma=" + gamma_str + "/Omega=" + Omega_str + "/g=" + g_str
        if not os.path.exists(file + subfile):
            os.makedirs(file + subfile)

    np.savetxt(file + "/output_data.txt", np.array(output_data), delimiter=',')