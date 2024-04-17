import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *
from numpy.fft import fft, ifft

def gaussian(x, A, x_0, sigma):
    return A * np.exp(- (x - x_0) ** 2 / (2 * sigma ** 2))


def resample(old_arrays, old_x, N_resample):
    new_arrays = []
    for i in range(len(old_arrays)):
        [new_array_i, new_x] = signal.resample(old_arrays[i], int(len(old_x) * (N_resample)), old_x)
        new_arrays.append(new_array_i)
    return new_arrays, new_x


if __name__ == '__main__':

    eq = 'pdnlS_nospace'
    t_rate = 1
    dt = 0.1
    T = 1000

    [alpha, beta, gamma_0, mu, nu] = [1, 1, 0.2, 0.1, 0.2]
    parameters = [alpha, beta, gamma_0, mu, nu]
    # Definición de la grilla
    [tmin, tmax, dt] = [0, T, dt]
    t_grid = np.arange(tmin, tmax + dt, dt)
    T = tmax
    Nt = t_grid.shape[0]

    # Initial Conditions Pattern
    U_1_init = 1 * (np.random.rand(1) + 1j * np.random.rand(1))
    #U_2_init = 0.01 * np.random.rand(1)
    operators = [0]

    # Empaquetamiento de parametros, campos y derivadas para integración
    fields_init = [U_1_init]
    grids = [t_grid, 0, 0]

    # Midiendo tiempo inicial
    now = datetime.datetime.now()
    print('Hora de Inicio: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
    time_init = time.time()

    final_fields, fields_history, time_grid = RK4_FD(eq, fields_init, parameters, grids, dt, Nt, operators, t_rate) #INTEGRACION EN EL TIEMPO

    now = datetime.datetime.now()
    print('Hora de Término: ' + str(now.hour) + ':' + str(now.minute) + ':' + str(now.second))
    time_fin = time.time()
    print(str(time_fin - time_init) + ' seg')

    # Reobteniendo campos
    U1_light = np.array(fields_history)[:, 0]
    #U2_light = np.array(fields_history)[:, 1]
    t_light = time_grid

    plt.plot(t_light, np.real(U1_light), color="b")
    plt.plot(t_light, np.imag(U1_light), color="r")
    plt.show()
    plt.close()
