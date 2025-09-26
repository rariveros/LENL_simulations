from fun import *

"""
===========================================================
  Objetivo:
 ---------
 Este código realiza simulaciones de osciladores de Duffing
 acoplados en redes de Erdős–Rényi. 

 - Genera y simula la dinámica de la red.
 - Guarda datos menos densos para graficar (U, V, fases).
 - Guarda datos más densos de la última parte de la simulación,
   destinados a calcular el espectro de Lyapunov.

 Dependencias:
 -------------
 Requiere el archivo `fun.py` con las funciones:
   - erdos_renyi_graph
   - RK4_FD
   - equations_FD
   - Der
===========================================================
"""

if __name__ == '__main__':

    ###########################################################
    # CONFIGURACIÓN GENERAL (cambiar rutas acorde a carpetas ordenadas)
    ###########################################################
    project_name = '/network_chimeras/FIG03/erdos_renyi'
    disc = 'D:/'
    route = 'mnustes_science/simulation_data/FD'

    # Parámetros de red a simular
    mean_degrees = [18, 22, 26]
    # Nota: para mean_degree >= 25 usar tmax=15000, para valores menores tmax=10000
    samples = [0]

    ###########################################################
    # LOOP SOBRE GRADOS PROMEDIO Y MUESTRAS
    ###########################################################
    for mean_degree in mean_degrees:
        print("########### MEAN DEGREE = " + str(mean_degree) + " ###########")

        for sample in samples:
            print("########### Sample = " + str(sample) + " ###########")
            now = datetime.datetime.now()
            print('Hora de Inicio: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
            time_init = time.time()

            ###########################################################
            # PARÁMETROS DE LA RED
            ###########################################################
            n = 501
            p = mean_degree / (n - 1)
            graph = erdos_renyi_graph(n, p)

            adj_matrix = nx.adjacency_matrix(graph).toarray()
            laplacian_matrix = nx.laplacian_matrix(graph).tocsc()
            L_dense = laplacian_matrix.toarray()

            ###########################################################
            # PARÁMETROS DEL OSCILADOR
            ###########################################################
            eq = 'duffing'  # Tipo de oscilador
            alpha = 0.4      # Coeficiente no lineal
            mu = 0.1         # Disipación
            gamma = 2.90     # Amplitud del forzamiento
            k = 0.015        # Acoplamiento
            w = 0.7          # Frecuencia del forzamiento
            t_rate = 1       # Frecuencia de guardado en simulación

            ###########################################################
            # PREPARACIÓN DE LA SIMULACIÓN
            ###########################################################
            N_nodes = n
            [tmin, tmax, dt] = [0, 5000, 0.05]   # Escala temporal
            t_grid = np.arange(tmin, tmax + dt, dt)   # Malla temporal

            [xmin, xmax, dx] = [0, N_nodes, 1]   # Escala espacial
            x_grid = np.arange(xmin, xmax, dx)   # Malla espacial

            Nt = t_grid.shape[0]
            Nx = x_grid.shape[0]

            # Condiciones iniciales
            U_init = 1.0 * np.ones(Nx)
            initial_quimera = 18
            arg_chimera = [initial_quimera]      # Índice inicial de quimera
            for i in arg_chimera:
                U_init[i] = 2.0

            U_init = U_init + 0.0 * (np.random.rand(Nx) - 0.5)
            V_init = 0.0 * np.random.rand(Nx)

            operators = [laplacian_matrix]
            fields_init = [U_init, V_init]
            grids = [t_grid, x_grid, 0]
            parameters_np = np.array([alpha, mu, gamma, k, w])

            ###########################################################
            # EJECUCIÓN DE LA SIMULACIÓN Y EMPAQUETADURA DE DATOS
            ###########################################################
            final_fields, fields_history, time_grid = RK4_FD(eq, fields_init, parameters_np, grids, dt, Nt, operators, t_rate)

            U = np.array(fields_history)[:, 0]
            V = np.array(fields_history)[:, 1]
            phase = np.arctan2(V, U)

            # Submuestreo (lightness)
            lightness = 1
            U_light = U[0::lightness]
            V_light = V[0::lightness]
            phase_light_wraped = phase[0::lightness]
            phase_light = np.unwrap(phase_light_wraped, axis=0)[0::lightness]
            t_light = np.array(time_grid[0::lightness])

            del fields_history
            gc.collect()

            ###########################################################
            # GUARDADO DE DATOS
            ###########################################################
            file = disc + route + project_name
            k_str = f"{k:.{4}f}"
            mean_degree_str = f"{mean_degree:.{2}f}"
            sample_str = f"{sample:.{2}f}"
            subfile = "/k=" + k_str + "/mean_degree=" + mean_degree_str

            if not os.path.exists(file + subfile):
                os.makedirs(file + subfile)

            # Datos menos densos (para graficar)
            np.savetxt(file + subfile + '/U.txt', U_light[::10], delimiter=',')
            np.savetxt(file + subfile + '/V.txt', V_light[::10], delimiter=',')
            np.savetxt(file + subfile + '/T.txt', t_light[::10], delimiter=',')
            np.savetxt(file + subfile + '/X.txt', x_grid, delimiter=',')
            np.savetxt(file + subfile + '/params.txt', parameters_np, delimiter=',')
            np.savetxt(file + subfile + '/L.txt', L_dense, delimiter=',')

            # Datos densos (para espectro de Lyapunov)
            N_init = int(0.9 * Nt)
            np.savetxt(file + subfile + '/U_lyap.txt', U_light[N_init:, :], delimiter=',')
            np.savetxt(file + subfile + '/V_lyap.txt', V_light[N_init:, :], delimiter=',')
            np.savetxt(file + subfile + '/T_lyap.txt', t_light[N_init:], delimiter=',')

            ###########################################################
            # LOG DE TIEMPOS
            ###########################################################
            now = datetime.datetime.now()
            print('Hora de Término: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
            time_fin = time.time()
            print(str(time_fin - time_init) + ' seg')

