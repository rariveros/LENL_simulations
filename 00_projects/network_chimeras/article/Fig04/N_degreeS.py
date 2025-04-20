# Importamos las librerías necesarias:
import numpy as np                       # Para cálculos numéricos y trabajar con arrays.
import matplotlib.pyplot as plt          # Para hacer gráficos.
import os                                # Para manejar carpetas y archivos del sistema.
import tkinter as tk                     # Para mostrar una ventanita.
from tkinter import filedialog           # Para permitir al usuario elegir una carpeta.
from back_process import *

# Solo se ejecuta si corres este script directamente (no si lo importas):
if __name__ == '__main__':

    frequencies = []  # Esta lista no se usa, pero está declarada...

    # Ruta inicial donde el usuario va a empezar a buscar
    disc = "D:/"  # Letra del disco (raro, pero probablemente es una unidad externa)
    initial_dir_data = str(disc) + 'mnustes_science/simulation_data/FD'

    # Mostramos una ventanita para que el usuario elija una carpeta
    root = tk.Tk()
    root.withdraw()  # Oculta la ventana principal de Tkinter
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')    #SE DEBE ELEGIR LA CARPETA EN DONDE ESTAN LOS DIFERENTES KAPPA

    save_directory = working_directory

    # Buscamos todas las subcarpetas dentro de la carpeta elegida
    directories = [name for name in os.listdir(working_directory) if os.path.isdir(os.path.join(working_directory, name))]

    # Preparamos una figura (gráfico) y un eje para dibujar
    fig, ax = plt.subplots(1, 1, figsize=(4, 3))  # 1 fila, 1 columna, tamaño de 5x3 pulgadas

    # Recorremos cada subcarpeta (cada una representa un valor distinto de kappa)
    for directory_01 in directories:
        dir_01 = working_directory + "/" + directory_01

        # Leemos el archivo 'data.txt' con los datos generados por el código anterior
        data = np.loadtxt(dir_01 + '/data.txt', delimiter=',')

        MEAN_DEGREES = data[0]   # Vector de grados medios de la red
        N_QUIM = data[1]         # Porcentaje promedio de nodos en estado quimera
        N_QUIM_std = data[2]     # Desviación estándar (qué tanto varía ese porcentaje)

        # Dibujamos los datos con barras de error (para mostrar la variabilidad)
        ax.errorbar(
            MEAN_DEGREES,        # Eje X: grado medio
            N_QUIM,              # Eje Y: fracción de quimeras
            yerr=N_QUIM_std,     # Barras de error en Y
            marker='o',          # Puntos como circulitos
            ls='',               # Sin línea entre puntos
            ecolor="k",          # Color negro para las líneas de error
            mec='black',         # Contorno negro para los marcadores
            label="$\kappa=" + directory_01.split('=')[-1] + "$"  # Leyenda con el valor de kappa sacado del nombre de la carpeta
        )

    # Ajustamos los límites de los ejes
    ax.set_xlim(13, 33)
    ax.set_ylim(-0.05, 1.05)

    # Definimos los valores del eje X manualmente
    ax.set_xticks([14, 16, 18, 20, 22, 24, 26, 28, 30, 32])

    # Personalizamos los "ticks" (las marcas en los ejes)
    ax.tick_params(axis="y", direction="in", labelsize=12, left=True, right=True, labelleft=True, labelright=False)
    ax.tick_params(axis="x", direction="in", labelsize=12, top=True, bottom=True, labeltop=False, labelbottom=True)

    # Etiquetas de los ejes (con formato matemático)
    ax.set_xlabel("$\\langle k \\rangle$", fontsize=20)      # Grado promedio
    ax.set_ylabel("$\\frac{N_c}{N}$", fontsize=20)           # Fracción de nodos en estado quimera

    # Textos extra en el gráfico (fijos)
    ax.text(13.5, 0.92, "$\kappa = 1.5 \\times 10^{-2}$", fontsize=12)  # Valor de kappa (puede que esté fijo, no automático)
    ax.text(13.5, 0.8, "$N = 501$", fontsize=12)                        # Tamaño de la red

    # Línea horizontal punteada en y = 0 y y = 1 como referencia
    ax.hlines([0, 1], 10, 35, linestyle="--", color="black")

    # Mostramos la leyenda (etiquetas de cada serie de datos)
    plt.legend()

    # Ajusta los márgenes automáticamente para que no se corte nada
    plt.tight_layout()

    # Guardamos el gráfico como imagen PNG
    plt.savefig('Fig03.png', dpi=200)
