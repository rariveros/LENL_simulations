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
            print(sample_folder_path)
            # Cargar archivos necesarios
            X_sorted_path = os.path.join(sample_folder_path, 'X_sorted.txt')
            chaos_path = os.path.join(sample_folder_path, 'chaos_no_chaos.txt')
            Adj_matrix_path = os.path.join(sample_folder_path, 'Adj_matrix.txt')

            # Cargar datos
            X_sorted = np.loadtxt(X_sorted_path).astype(np.int64)  # array de enteros
            chaos_no_chaos = np.loadtxt(chaos_path)

            # Reordenar
            chaos_ordered = chaos_no_chaos#[X_sorted]
            #print(X_sorted)
            CHAOS.append(chaos_ordered)
    Adj_matrix = np.loadtxt(Adj_matrix_path, delimiter=",")
    G = nx.from_numpy_array(Adj_matrix)
    degree_centrality = nx.closeness_centrality(G)
    CHAOS = np.array(CHAOS)
    CHAOS_mean = np.mean(CHAOS, axis=0)
    CHAOS_std = np.std(CHAOS, axis=0)
    #print(degree_centrality)
    #print(CHAOS_mean)
    plt.scatter(degree_centrality.values(), CHAOS_mean)
    plt.xlabel("degree")
    plt.ylabel("chaos")
    #plt.xlim(0, np.amax(degree_centrality.values()))
    plt.ylim(0, np.amax(CHAOS_mean))
    #plt.scatter(np.arange(0, 501), np.sort(CHAOS_mean))
    #plt.errorbar(np.arange(0,501), CHAOS_mean, CHAOS_std, fmt="o")
    plt.show()
