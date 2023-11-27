from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':

    # Definiendo parámetros
    project_name = '/PDNLS_cos_num'
    disc = 'C:/'                                        # DISCO DE TRABAJO
    route = 'mnustes_science/simulation_data/FD'        # CARPETA DE TRABAJO
    eq = 'transcritical'                                        # ECUACION
    t_rate = 2                                        # CADA CUANTAS ITERACIONES GUARDA
    dt = 0.0001
    T = 200

    alphas = np.arange(-0.1, 0.1, 0.025)
    finals = []
    for alpha in alphas:
        finals_i = []
        # Midiendo tiempo inicial
        now = datetime.datetime.now()
        print('Hora de Inicio: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
        time_init = time.time()
        for i in range(0, 3):
            #alpha = -0.02  #5.721
            beta = 0.1
            gamma = 0

            # Definición de la grilla
            [tmin, tmax, dt] = [0, T, dt]
            t_grid = np.arange(tmin, tmax + dt, dt)
            T = tmax
            Nt = t_grid.shape[0]

            # Initial Conditions Pattern
            U_1_init = 0.1 * np.array([1])

            # Empaquetamiento de parametros, campos y derivadas para integración
            fields_init = [U_1_init]
            grids = [t_grid, 0, 0]
            phi = np.pi

            parameters = [alpha, beta, gamma]

            final_fields, fields_history, time_grid = RK4_ode(eq, fields_init, t_grid, parameters, dt, Nt)

            # Reobteniendo campos
            U1_light = np.array(fields_history)[:, 0]
            t_light = time_grid
            finals_i.append(U1_light[-1])
            #plt.plot(t_light, U1_light)
            #plt.show()
        now = datetime.datetime.now()
        print('Hora de Término: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
        time_fin = time.time()
        print(str(time_fin - time_init) + ' seg')
        average = np.mean(np.array(finals_i))
        finals.append(average)
    print(np.array(finals))
    plt.scatter(alphas, np.array(finals))
    plt.show()

    # Definiendo variables finales
    #modulo_light_1 = np.absolute(U_complex)
    #arg_light_1 = np.angle(U_complex)
    #arg_light_1 = (2 * np.pi + arg_light_1) * (arg_light_1 < 0) + arg_light_1 * (arg_light_1 > 0)
    #analytical_signal_1 = hilbert(U1_light[-1, :])
    #amplitude_envelope_1 = np.abs(analytical_signal_1)

    # Guardando datos
    #file = disc + route + project_name
    #subfile = pdnlS_bigauss_name([alpha, beta, mu, nu, sigma, gamma_0, dist], "ABMNSGD")
    #parameters_np = np.array([alpha, beta, gamma_0, mu, nu, sigma])

    #if not os.path.exists(file + subfile):
    #    os.makedirs(file + subfile)
    #np.savetxt(file + subfile + '/field_real.txt', U1_light, delimiter=',')
    #np.savetxt(file + subfile + '/field_img.txt', U2_light, delimiter=',')
    #np.savetxt(file + subfile + '/parameters.txt', parameters_np, delimiter=',')
    #np.savetxt(file + subfile + '/final_envelope.txt', amplitude_envelope_1, delimiter=',')
    #np.savetxt(file + subfile + '/X.txt', x_grid, delimiter=',')
    #np.savetxt(file + subfile + '/T.txt', t_light, delimiter=',')