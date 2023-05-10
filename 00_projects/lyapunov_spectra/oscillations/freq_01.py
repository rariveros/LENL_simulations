import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    disco = 'E:/'
    initial_dir_data = str(disco) + 'Users/mnustes_science/PT_fluids/mnustes_science/simulation_data/simulation_data'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    Z_r = np.loadtxt(directory + '/Z_r_points.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/Z_i_points.txt', delimiter=',')
    X = np.loadtxt(directory + '/x_points.txt', delimiter=',')
    T = np.loadtxt(directory + '/t_points.txt', delimiter=',')

    Nt = len(T)
    dt = np.abs(T[1] - T[0])

    yf_i = fft(Z_r)
    yf_i = 2.0 / Nt * np.abs(yf_i[0:Nt // 2])
    xf_i = fftfreq(Nt, dt)[:Nt // 2]

    plt.plot(xf_i[1:], yf_i[1:])
    plt.show()
    plt.close()