import networkx as nx
import matplotlib.pyplot as plt
import numpy as np


if __name__ == '__main__':
    # Crear red BA
    G = nx.LFR_benchmark_graph(
        n=1000,
        tau1=3,  # Exponente de distribución de grados
        tau2=2,  # Exponente de distribución del tamaño de comunidades
        mu=0.3,  # Mezcla: 0 = comunidades claras, 1 = sin estructura
        average_degree=8,
        min_community=10,
        seed=42
    )

    # Obtener grados
    degrees = [d for n, d in G.degree()]

    # Histograma en escala log-log
    plt.figure(figsize=(6,4))
    hist, bins = np.histogram(degrees, bins=range(min(degrees), max(degrees)+1))
    plt.scatter(bins[:-1], hist)
    plt.xlabel('Grado k')
    plt.ylabel('P(k)')
    plt.title('Distribución de grados en red BA')
    plt.grid(True)
    plt.show()