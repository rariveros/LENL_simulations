from fun import *
"""
===========================================================
 Objetivo:
 ---------
 Este código simula osciladores de Duffing en una red 
 Erdős–Rényi fija, explorando la dinámica al variar las 
 condiciones iniciales.

 - Genera una sola red al inicio de la simulación.
 - Guarda la matriz de adyacencia de la red en el nivel 
   superior de carpetas.
 - Corre simulaciones para distintos valores de "sample"
   que representan distintas condiciones iniciales.
 - Calcula promedios temporales de módulo y fase por nodo
   usando el último 10% del tiempo de simulación.
 - Ordena nodos según el módulo promedio y guarda el vector
   de índices de ordenamiento (args).
 - Guarda resultados en orden natural (índice en x).
 - Genera figuras espacio–temporales antes y después del
   reordenamiento.

 Nota:
 -----
 Si en vez de generar una red quieres **cargar una matriz
 de adyacencia** ya existente, descomenta el bloque marcado
 como "OPCIÓN: CARGAR MATRIZ DE ADJACENCIA".
===========================================================
"""

if __name__ == '__main__':

    ###########################################################
    # CONFIGURACIÓN GENERAL
    ###########################################################
    project_name = '/network_chimeras/FIGxx/erdos_renyi_CI'
    disc = 'D:/'
    route = 'mnustes_science/simulation_data/FD'

    Ks = [0.016]          # valor fijo de k
    samples = np.arange(0, 20)   # condiciones iniciales distintas

    # Parámetros generales de la red
    mean_degree = 20
    n = 501

    ###########################################################
    # CREAR O CARGAR LA RED UNA SOLA VEZ
    ###########################################################
    p = mean_degree / (n - 1)
    graph = erdos_renyi_graph(n, p)
    adj_matrix = nx.adjacency_matrix(graph).toarray()
    degrees = np.sum(adj_matrix, axis=1)
    D = np.diag(degrees)

    # Matriz laplaciana
    L_dense = D - adj_matrix
    laplacian_matrix = L_dense

    # Guardar matriz de adyacencia una vez, al nivel superior
    file = disc + route + project_name
    if not os.path.exists(file):
        os.makedirs(file)
    np.savetxt(file + '/Adj_matrix.txt', adj_matrix, delimiter=',')

    # =========================================================
    # OPCIÓN: CARGAR MATRIZ DE ADJACENCIA
    # adj_matrix = np.loadtxt("ruta/a/Adj_matrix.txt", delimiter=',')
    # degrees = np.sum(adj_matrix, axis=1)
    # D = np.diag(degrees)
    # L_dense = D - adj_matrix
    # laplacian_matrix = L_dense
    # =========================================================

    ###########################################################
    # LOOP SOBRE k Y CONDICIONES INICIALES
    ###########################################################
    for k in Ks:
        print("########### k = " + str(k) + " ###########")

        for sample in samples:
            print("########### Sample = " + str(sample) + " ###########")
            now = datetime.datetime.now()
            print('Hora de Inicio: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
            time_init = time.time()

            ###########################################################
            # PARÁMETROS DEL OSCILADOR
            ###########################################################
            alpha = 0.4      # Coeficiente no lineal
            mu = 0.1         # Disipación
            gamma = 2.90     # Fuerza de forzamiento
            w = 0.7          # Frecuencia de forzamiento
            eq = 'duffing'
            t_rate = 5       # Frecuencia de guardado

            ###########################################################
            # PREPARACIÓN DE LA SIMULACIÓN
            ###########################################################
            N_nodes = n
            [tmin, tmax, dt] = [0, 100, 0.05]
            t_grid = np.arange(tmin, tmax + dt, dt)    # Malla temporal

            [xmin, xmax, dx] = [0, N_nodes, 1]
            x_grid = np.arange(xmin, xmax, dx)         # Malla espacial

            Nt = t_grid.shape[0]
            Nx = x_grid.shape[0]

            # Condiciones iniciales distintas por cada sample
            U_init = 1.0 * np.ones(Nx)
            initial_quimera = np.random.randint(0, N_nodes)
            U_init[initial_quimera] = 2.0
            U_init = U_init + 0.1 * (np.random.rand(Nx) - 0.5)
            V_init = 0.1 * (np.random.rand(Nx) - 0.5)

            operators = [laplacian_matrix]
            fields_init = [U_init, V_init]
            grids = [t_grid, x_grid, 0]
            parameters_np = np.array([alpha, mu, gamma, k, w])

            ###########################################################
            # EJECUCIÓN DE LA SIMULACIÓN
            ###########################################################
            final_fields, fields_history, time_grid = RK4_FD(eq, fields_init, parameters_np, grids, dt, Nt, operators, t_rate )
            time_grid = np.array(time_grid)

            ###########################################################
            # PROCESAMIENTO DE RESULTADOS
            ###########################################################
            # Usar el último 10% del tiempo
            N_init = int(0.9 * len(time_grid))
            U_light = np.array(fields_history)[N_init:, 0]
            V_light = np.array(fields_history)[N_init:, 1]
            t_light = np.array(time_grid[N_init:])

            del fields_history
            gc.collect()

            # Calcular módulo y fase
            phase_light_wraped = np.arctan2(V_light, U_light)
            phase_light = np.unwrap(phase_light_wraped, axis=0)
            module = np.sqrt(U_light ** 2 + V_light ** 2)

            # Promedios temporales
            average_phase = integrate.simpson(phase_light_wraped, t_light, axis=0) / (t_light[-1] - t_light[0])
            average_module = integrate.simpson(module, t_light, axis=0) / (t_light[-1] - t_light[0])

            # Ordenar nodos por módulo promedio
            args = np.argsort(average_module)

            ###########################################################
            # GUARDADO DE RESULTADOS
            ###########################################################
            k_str = f"{k:.{4}f}"
            mean_degree_str = f"{mean_degree:.{2}f}"
            sample_str = f"{sample:.{2}f}"
            subfile = "/k=" + k_str + "/mean_degree=" + mean_degree_str + "/sample=" + sample_str

            if not os.path.exists(file + subfile):
                os.makedirs(file + subfile)

            # Guardar datos en orden natural
            np.savetxt(file + subfile + '/module_mean.txt', average_module, delimiter=',')
            np.savetxt(file + subfile + '/phase_mean.txt', average_phase, delimiter=',')
            np.savetxt(file + subfile + '/X.txt', x_grid, delimiter=',')

            # Guardar args para reordenar externamente
            np.savetxt(file + subfile + '/args_module.txt', args, delimiter=',', fmt='%d')

            # Guardar parámetros
            np.savetxt(file + subfile + '/params.txt', parameters_np, delimiter=',')

            ###########################################################
            # FIGURA ESPACIO-TEMPORAL
            ###########################################################
            # Seleccionar última ventana pequeña para graficar
            t_init_fig, t_final_fig = 9950, 10000
            i_0 = np.argmin(np.abs(t_light - t_init_fig))
            i_f = np.argmin(np.abs(t_light - t_final_fig))

            fig, (ax01, ax02) = plt.subplots(2, 1, figsize=(7, 6))

            # Orden natural
            cax_01 = ax01.pcolormesh(x_grid, t_light[i_0:i_f], module[i_0:i_f, :], cmap="turbo", shading='auto')
            fig.colorbar(cax_01, ax=ax01, label='$r(t)$')
            ax01.set_xlabel("$i$", fontsize=20)
            ax01.set_ylabel("$t$", fontsize=20)

            # Ordenado por módulo promedio
            cax_02 = ax02.pcolormesh(x_grid, t_light[i_0:i_f], module[i_0:i_f, args], cmap="turbo", shading='auto')
            fig.colorbar(cax_02, ax=ax02, label='$r(t)$')
            ax02.set_xlabel("$i$", fontsize=20)
            ax02.set_ylabel("$t$", fontsize=20)

            plt.tight_layout()
            plt.savefig(file + subfile + '/spatiotemporal.png', dpi=200)
            plt.close(fig)

            ###########################################################
            # LOG DE TIEMPOS
            ###########################################################
            now = datetime.datetime.now()
            print('Hora de Término: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
            time_fin = time.time()
            print(str(time_fin - time_init) + ' seg')