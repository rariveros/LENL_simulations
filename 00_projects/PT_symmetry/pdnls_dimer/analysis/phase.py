import numpy as np

from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    disc = "aD:"
    initial_dir_data = str(disc) + 'mnustes_science/simulation_data/FD'
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    dir_01 = working_directory
    params = np.loadtxt(dir_01 + '/parameters.txt', delimiter=',')
    T = np.loadtxt(dir_01 + '/T.txt', delimiter=',')
    U1 = np.loadtxt(dir_01 + '/U1.txt', delimiter=',', dtype=np.complex128)
    V1 = np.loadtxt(dir_01 + '/V1.txt', delimiter=',', dtype=np.complex128)
    save_directory = working_directory
    Nt = len(T)
    t0 = int(0.5 * Nt)

    T = T[t0:] - T[t0]
    U1 = U1[t0:]
    V1 = V1[t0:]
    Nt = len(T)
    dt = T[1] - T[0]

    arg_u = np.arctan2(np.real(U1), np.imag(U1))
    arg_v = np.arctan2(np.real(V1), np.imag(V1))
    plt.scatter(np.unwrap(arg_u) + np.pi, np.unwrap(arg_v) + np.pi, color="k")
    plt.show()

