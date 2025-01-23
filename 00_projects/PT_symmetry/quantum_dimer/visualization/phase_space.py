from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    working_directory = "D:/mnustes_science/simulation_data/FD/ladder_operators/bifurcation_steady/Delta=0.1000/gamma=0.1000"
    directories = [name for name in os.listdir(working_directory) if os.path.isdir(os.path.join(working_directory, name))]

    for directory in directories:
        U1 = np.loadtxt(working_directory + "/" + directory + '/k=0.2600/U1.txt', delimiter=',', dtype=np.complex128)
        V1 = np.loadtxt(working_directory + "/" + directory + '/k=0.2600/V1.txt', delimiter=',', dtype=np.complex128)
        T = np.loadtxt(working_directory + "/" + directory + '/k=0.2600/T.txt', delimiter=',')
        Nt = len(T)
        t0 = int(0.5 * Nt)
        T = T[t0:] - T[t0]
        U1 = U1[t0:]
        V1 = V1[t0:]
        Nt = len(T)
        dt = T[1] - T[0]
        plt.plot(np.imag(U1), np.real(V1), color="k")
        #plt.plot(-np.imag(U1),-np.real(V1), color="k")
        #plt.plot(np.imag(V1 * 1j), np.real(U1 * 1j), color="k")
        #plt.plot(-np.imag(V1 * 1j), -np.real(U1 * 1j), color="k")
    plt.show()
