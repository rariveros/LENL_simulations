import os
import re
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import gaussian_kde
import tkinter as tk
from tkinter import filedialog
from scipy.optimize import curve_fit


# =======================
# FUNCIONES AUXILIARES PARA PLOTEO
# =======================

def plot_density(average_module, k, mean_degree, sample, output_folder):
    """Genera y guarda el histograma con densidad de los valores average_module."""
    density = gaussian_kde(average_module)
    xs = np.linspace(min(average_module), max(average_module), 200)
    density_vals = density(xs)

    plt.figure(figsize=(8, 6))
    plt.hist(average_module, bins=50, color='skyblue', alpha=0.6, density=True)
    plt.plot(xs, density_vals, color='darkblue', linewidth=2)
    plt.title(f"Histograma y densidad\nk={k}, mean_degree={mean_degree}, sample={sample}")
    plt.xlabel("average_module")
    plt.ylabel("Densidad")
    plt.grid(True)
    plt.tight_layout()

    os.makedirs(output_folder, exist_ok=True)
    filename = f"k={k}_mean_degree={mean_degree}_sample={sample}.png"
    plt.savefig(os.path.join(output_folder, filename))
    plt.close()


def cumulate_function(average_module):
    kde = gaussian_kde(average_module)
    x_vals = np.linspace(min(average_module), max(average_module), 1000)
    f_x = kde(x_vals)

    dx = x_vals[1] - x_vals[0]
    F_x = np.cumsum(f_x) * dx
    F_x /= F_x[-1]  # Normalizar a [0, 1]

    return x_vals, F_x


def plot_cdf(x_vals, F_x, sample, output_folder):
    """Genera y guarda la función de distribución acumulada estimada de average_module."""

    plt.figure(figsize=(8, 6))
    plt.plot(x_vals, F_x, color='green', label='CDF estimada')
    plt.title(f'CDF acumulada\nk={k}, mean_degree={mean_degree}, sample={sample}')
    plt.xlabel("average_module")
    plt.ylabel("Probabilidad acumulada")
    plt.grid(True)
    plt.legend()

    os.makedirs(output_folder, exist_ok=True)
    filename = f"k={k}_mean_degree={mean_degree}_sample={sample}.png"
    plt.savefig(os.path.join(output_folder, filename))
    plt.close()


# Función modelo definida externamente
def model_func(x, a, xc, delta_c, delta_xs, delta_s):
    xs = xc + delta_xs  # impone la restricción xc < xs
    return 0.5 * (a * np.tanh((x - xc) / delta_c) + (1 - a) * np.tanh((x - xs) / delta_s)) + 0.5


# Función de ajuste
def ajustar_curva(x_vals, F_x):
    x_vals = np.array(x_vals)
    F_x = np.array(F_x)
    x_max = x_vals[-1]
    x_min = x_vals[0]

    # # Límites para los parámetros: a, xc, delta_c, delta_xs, delta_s
    # bounds = (
    #     [0,     x_min,     1e-6,     0,  1e-6],     # Inferiores
    #     [1, x_max, np.inf, x_max-x_min,  np.inf]             # Superiores
    # )

    # # Estimación inicial razonable: [a, xc, delta_c, delta_xs, delta_s]
    # initial_guess = [0.5, x_min+0.1, (x_max-x_min)/4, (x_max-x_min)*0.9, (x_max-x_min)/4]

    # Límites para los parámetros: a, xc, delta_c, delta_xs, delta_s
    bounds = (
        [0, x_min, 1e-6, 0.2, 0.01],  # Inferiores
        [1, 1.5, np.inf, 1, 2]  # Superiores
    )

    # Estimación inicial razonable: [a, xc, delta_c, delta_xs, delta_s]
    initial_guess = [0.5, x_min + 0.01, (x_max - x_min) / 4, 0.6, (x_max - x_min) / 4]

    # Ajustar
    popt, _ = curve_fit(model_func, x_vals, F_x, p0=initial_guess, bounds=bounds)

    return popt  # [a, xc, delta_c, delta_xs, delta_s]


# =======================
# EJECUCIÓN PRINCIPAL
# =======================
if __name__ == '__main__':

    # Configurar ruta inicial
    disc = "/Users/danieltoro/"
    initial_dir_data = os.path.join(disc, 'Data/erdos_renyi')

    # Diálogo para seleccionar carpeta raíz
    root = tk.Tk()
    root.withdraw()
    root_directory = filedialog.askdirectory(
        parent=root,
        initialdir=initial_dir_data,
        title='Selecciona la carpeta raíz de todas las simulaciones'
    )

    results = []
    fig, ax = plt.subplots(1, 1, figsize=(5, 3))

    # Iterar sobre subdirectorios
    for k_folder_name in os.listdir(root_directory):
        k_folder_path = os.path.join(root_directory, k_folder_name)
        if not os.path.isdir(k_folder_path):
            continue

        print(f"\n[INFO] Procesando carpeta: {k_folder_name}")
        Nc_means = []
        for mean_degree_folder_name in os.listdir(k_folder_path):
            mean_degree_folder_path = os.path.join(k_folder_path, mean_degree_folder_name)
            if not os.path.isdir(mean_degree_folder_path):
                continue
            Nc = []
            for sample_folder_name in os.listdir(mean_degree_folder_path):
                sample_folder_path = os.path.join(mean_degree_folder_path, sample_folder_name)
                if not os.path.isdir(sample_folder_path):
                    continue

                # Cargar archivos necesarios
                module_mean_path = os.path.join(sample_folder_path, 'module_mean.txt')
                phase_mean_path = os.path.join(sample_folder_path, 'phase_mean.txt')

                match = re.search(r'k=([0-9.]+)\\mean_degree=([0-9.]+)\\sample=([0-9.]+)', sample_folder_path)
                if match:
                    k = float(match.group(1))
                    mean_degree = float(match.group(2))
                    sample = float(match.group(3))
                else:
                    print(f"[WARN] No se pudo extraer metadata de la ruta: {sample_folder_path}")
                    continue

                if not (os.path.exists(module_mean_path) and os.path.exists(phase_mean_path)):
                    print(f"[WARN] Archivos no encontrados en {sample_folder_path}. Saltando...")
                    continue

                # Leer datos
                average_module = np.loadtxt(module_mean_path, delimiter=',')
                x_sorted = np.loadtxt(os.path.join(sample_folder_path, 'X_sorted.txt'), delimiter=',')

                # ajuste de la distribución a una función acumulada
                x_vals, F_x = cumulate_function(average_module)

                # ajustamos esta curva a una función
                popt = ajustar_curva(x_vals, F_x)
                a, xc, delta_c, delta_xs, delta_s = popt

                F_x_fit = model_func(x_vals, a, xc, delta_c, delta_xs, delta_s)


                xs = xc + delta_xs

                #plt.plot(x_vals, F_x, color="k")
                #plt.plot(x_vals, F_x_fit, color="r", ls="--")
                #plt.savefig(os.path.join(sample_folder_path, 'cdf_fit.png'), dpi=200)
                #plt.close()
                module_cut = (xc + xs) / 2
                chaos_no_chaos = []
                for i in range(len(average_module)):
                    if average_module[i] > module_cut:
                        chaos_no_chaos.append(0)
                    else:
                        chaos_no_chaos.append(1)
                chaos_no_chaos = np.array(chaos_no_chaos)
                np.savetxt(os.path.join(sample_folder_path, 'chaos_no_chaos.txt'), chaos_no_chaos, delimiter=',')
                # cuimera size

                chim_size = np.sum(chaos_no_chaos) / len(chaos_no_chaos)#a

                # save the observation

                observation = [k, mean_degree, sample, a, xc, delta_c, xs, delta_s]
                results.append(observation)

                # Imprimir los parámetros ajustados
                print(f"a = {a:.4f}, xc = {xc:.4f}, delta_c = {delta_c:.4f}, xs = {xs:.4f}, delta_s = {delta_s:.4f}")
                # Graficar los datos originales y el ajuste
                x_fit = np.linspace(min(x_vals), max(x_vals), 500)
                y_fit = model_func(x_fit, *popt)
                Nc.append(chim_size)
            Nc_means = np.mean(np.array(Nc))
            Nc_std = np.std(np.array(Nc))
            ax.errorbar(k, Nc_means, Nc_std, marker='o', ls='', ecolor="k", mec='black', color="k")
    #ax.set_xlim(13, 33)
    ax.set_ylim(-0.05, 1.05)
    #ax.set_xticks([14, 16, 18, 20, 22, 24, 26, 28, 30, 32])
    ax.tick_params(axis="y", direction="in", labelsize=12, left=True, right=True, labelleft=True,
                   labelright=False)
    ax.tick_params(axis="x", direction="in", labelsize=12, top=True, bottom=True, labeltop=False,
                   labelbottom=True)
    ax.set_xlabel("$\\langle k \\rangle$", fontsize=20)
    ax.set_ylabel("$\\frac{N_c}{N}$", fontsize=20)
    ax.text(13.5, 0.92, "$\kappa = 1.5 \\times 10^{-2}$", fontsize=12)
    ax.text(13.5, 0.8, "$N = 501$", fontsize=12)
    plt.tight_layout()
    plt.savefig('Fig03.png', dpi=200)
    plt.show()
    plt.close()

    df = pd.DataFrame(results, columns=["k", "mean_degree", "sample", "a", "xc", "delta_c", "xs", "delta_s"])
    df.to_csv(os.path.join(root_directory, "resultados.csv"), index=False, sep=";")

    # DANIEL QUIERE: INDEX, QUIMERA NO QUIMERA, GRAFO