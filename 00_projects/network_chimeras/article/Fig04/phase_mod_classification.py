# Importamos librerías necesarias:
import numpy as np             # Para trabajar con arreglos y cálculos numéricos.
import os                      # Para manejar carpetas y archivos del sistema operativo.
import tkinter as tk           # Para hacer una ventanita de selección de carpeta.
from tkinter import filedialog # Para que el usuario pueda elegir una carpeta.

# Esto solo se ejecuta si el archivo se corre directamente (no si lo importas como módulo):
if __name__ == '__main__':

    # Ruta donde parte la selección de carpetas
    disc = "C:"  # Disco duro donde están las carpetas
    initial_dir_data = str(disc) + 'mnustes_science/simulation_data/FD'  # Ruta inicial de datos

    # Se crea una ventanita para que el usuario elija una carpeta de trabajo (debe ser por cada k de kappa de acople)
    root = tk.Tk()
    root.withdraw()  # Oculta la ventana principal de Tkinter
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    # Definimos la carpeta de guardado de los datos
    save_directory = working_directory

    # Dentro de esa carpeta buscamos todas las subcarpetas
    directories = [name for name in os.listdir(working_directory) if os.path.isdir(os.path.join(working_directory, name))]

    # Inicializamos tres listas para guardar los resultados
    MEAN_DEGREES = []  # Aquí se guarda el grado medio de la red
    N_QUIM = []        # Aquí se guarda el porcentaje de osciladores en estado quimera
    N_QUIM_std = []    # Aquí se guarda la desviación estándar del porcentaje de quimeras

    # Recorremos cada subcarpeta (cada una representa una simulación distinta)
    for directory_01 in directories:
        dir_01 = working_directory + "/" + directory_01
        directories_01 = [name for name in os.listdir(dir_01) if os.path.isdir(os.path.join(dir_01, name))]     # Buscamos más subcarpetas dentro de cada carpeta principal
        n_quim = []         # Para guardar el % de quimeras de cada simulación individual
        mean_degree = []    # Para guardar el grado medio de cada simulación

        for directory_02 in directories_01:
            dir_02 = working_directory + "/" + directory_01 + "/" + directory_02

            # Leemos dos archivos que tienen datos de la simulación
            output = np.loadtxt(dir_02 + '/output.txt', delimiter=',')           # output.txt tiene información general
            average_module = np.loadtxt(dir_02 + '/module_mean.txt', delimiter=',') # Este tiene la medida de sincronización de cada nodo

            N_nodes = len(average_module)  # Número de nodos en la red

            # Según un valor de output, se elige un umbral para decidir si un nodo está en estado quimera o no (ACA HAY QUE COLOCAR EL CLASIFICADOR SINCR/QUIM)
            if output[1] <= 25:
                power_threshold = 1.5
            else:
                power_threshold = 1.6

            # Separamos los nodos: los que están en quimera (por debajo del umbral) y los sincronizados
            arg_chimeras = average_module < power_threshold
            arg_sync = average_module >= power_threshold

            # Calculamos el porcentaje de nodos en estado quimera
            frac_chimera = len(average_module[arg_chimeras]) / N_nodes

            # Guardamos el porcentaje de quimeras y el grado medio
            n_quim.append(frac_chimera)
            mean_degree.append(output[1])

        # Guardamos los resultados para esa red (una por carpeta grande)
        MEAN_DEGREES.append(mean_degree[0])     # Se apendean los degrees
        N_QUIM.append(np.mean(n_quim))          # Promedio de % de quimeras
        N_QUIM_std.append(np.std(n_quim))       # Guardamos esa desviación estándar

    # Convertimos tutto a arreglos de numpy para guardar en un archivo
    MEAN_DEGREES = np.array(MEAN_DEGREES)
    N_QUIM = np.array(N_QUIM)
    N_QUIM_std = np.array(N_QUIM_std)

    # Guardamos los datos en un archivo de texto llamado "data.txt"
    np.savetxt(save_directory + '/data.txt', [MEAN_DEGREES, N_QUIM, N_QUIM_std], delimiter=',')
