import os
import re
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import gaussian_kde
import tkinter as tk
from tkinter import filedialog
from scipy.optimize import curve_fit
import networkx as nx

from back_process import *

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

    CHAOS = []

    # Iterar sobre subdirectorios
    Nc_means = []
    for mean_degree_folder_name in os.listdir(root_directory):
        mean_degree_folder_path = os.path.join(root_directory, mean_degree_folder_name)
        if not os.path.isdir(mean_degree_folder_path):
            continue
        Nc = []
        for sample_folder_name in os.listdir(mean_degree_folder_path):
            sample_folder_path = os.path.join(mean_degree_folder_path, sample_folder_name)
            if not os.path.isdir(sample_folder_path):
                continue

            # Cargar archivos necesarios
            X_sorted_path = os.path.join(sample_folder_path, 'args_module.txt')
            chaos_path = os.path.join(sample_folder_path, 'chaos_no_chaos.txt')
            Adj_matrix_path = os.path.join(r"C:\mnustes_science\simulation_data\FD\NW_chimeras\FIGxx\erdos_renyi_CI", 'Adj_matrix.txt')

            # Cargar datos
            X_sorted = np.loadtxt(X_sorted_path, delimiter=',', dtype="int")  # array de enteros
            chaos_no_chaos = np.loadtxt(chaos_path)  # array de floats (o str si quieres dtype=str)
            chaos_ordered = chaos_no_chaos # reorder chaos array the same ways
            CHAOS.append(chaos_ordered)
    Adj_matrix = np.loadtxt(Adj_matrix_path, delimiter=",")
    G = nx.from_numpy_array(Adj_matrix)
    degree_centrality = nx.degree_centrality(G)
    betweenness_centrality = nx.betweenness_centrality(G, normalized=True)
    eigenvector_centrality = nx.eigenvector_centrality_numpy(G)
    CHAOS = np.array(CHAOS)
    CHAOS_mean = np.mean(CHAOS, axis=0)
    CHAOS_std = np.std(CHAOS, axis=0)

    # =========================================================
    # === FIGURA CON 3 HISTOGRAMAS (UNA POR CENTRALIDAD) ===
    # =========================================================

    # ---- Degree ----
    x_deg = np.array(list(degree_centrality.values()))
    y = CHAOS_mean
    bins_deg = np.linspace(x_deg.min(), x_deg.max(), 20)
    bin_centers_deg = 0.5 * (bins_deg[:-1] + bins_deg[1:])
    bin_means_deg = [y[(x_deg >= bins_deg[i]) & (x_deg < bins_deg[i + 1])].mean() for i in range(len(bins_deg) - 1)]
    bin_stds_deg = [y[(x_deg >= bins_deg[i]) & (x_deg < bins_deg[i + 1])].std() for i in range(len(bins_deg) - 1)]

    # ---- Betweenness ----
    x_bet = np.array(list(betweenness_centrality.values()))
    bins_bet = np.linspace(x_bet.min(), x_bet.max(), 20)
    bin_centers_bet = 0.5 * (bins_bet[:-1] + bins_bet[1:])
    bin_means_bet = [y[(x_bet >= bins_bet[i]) & (x_bet < bins_bet[i + 1])].mean() for i in range(len(bins_bet) - 1)]
    bin_stds_bet = [y[(x_bet >= bins_bet[i]) & (x_bet < bins_bet[i + 1])].std() for i in range(len(bins_bet) - 1)]

    # ---- Eigenvector ----
    x_eig = np.array(list(eigenvector_centrality.values()))
    bins_eig = np.linspace(x_eig.min(), x_eig.max(), 20)
    bin_centers_eig = 0.5 * (bins_eig[:-1] + bins_eig[1:])
    bin_means_eig = [y[(x_eig >= bins_eig[i]) & (x_eig < bins_eig[i + 1])].mean() for i in range(len(bins_eig) - 1)]
    bin_stds_eig = [y[(x_eig >= bins_eig[i]) & (x_eig < bins_eig[i + 1])].std() for i in range(len(bins_eig) - 1)]

    # =========================================================
    # === PLOT (3 SUBFIGURAS HORIZONTALES) ===
    # =========================================================
    fig, (ax01, ax02, ax03) = plt.subplots(1, 3, figsize=(5.5, 2.5))

    # --- Degree ---
    ax01.bar(bin_centers_deg, bin_means_deg,
             width=(bins_deg[1] - bins_deg[0]) * 1,
             yerr=bin_stds_deg, capsize=2,
             color='cornflowerblue', edgecolor='k', ecolor='k', alpha=0.85)
    ax01.set_xlabel(r"\textrm{Degree}", fontsize=13)
    ax01.set_ylabel(r"$P_{\textrm{chaos}}$", fontsize=13)
    ax01.tick_params(axis="both", direction="in", labelsize=11, top=True, bottom=True, right=True, left=False)
    ax01.set_ylim(0, 1.1)

    # --- Betweenness ---
    ax02.bar(bin_centers_bet, bin_means_bet,
             width=(bins_bet[1] - bins_bet[0]) * 1,
             yerr=bin_stds_bet, capsize=2,
             color='indianred', edgecolor='k', ecolor='k', alpha=0.85)
    ax02.set_xlabel(r"\textrm{Betweenness}", fontsize=13)
    ax02.tick_params(axis="both", direction="in", labelsize=11, top=True, bottom=True, right=True, left=False, labelleft=False)
    ax02.set_ylim(0, 1.1)

    # --- Eigenvector ---
    ax03.bar(bin_centers_eig, bin_means_eig,
             width=(bins_eig[1] - bins_eig[0]) * 1,
             yerr=bin_stds_eig, capsize=2,
             color='mediumseagreen', edgecolor='k', ecolor='k', alpha=0.85)
    ax03.set_xlabel(r"\textrm{Eigenvector}", fontsize=13)
    ax03.tick_params(axis="both", direction="in", labelsize=11, top=True, bottom=True, right=True, left=False, labelleft=False)
    ax03.set_ylim(0, 1.1)

    # Ajustes de layout (horizontal)
    fig.subplots_adjust(wspace=0.1, hspace=0.0, left=0.13, right=0.98, top=0.93, bottom=0.18)
    plt.savefig("centrality_histograms_binned.png", dpi=300)
    # plt.show()



