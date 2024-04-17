import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    disco = 'E:/'
    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FD'
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    directories = os.listdir(working_directory)
    arg_max = []
    mod_max = []
    x_arg_max = []
    x_mod_max = []
    #fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    N_dir = len(directories)
    sigmas = np.arange(4, 25, 1)
    for i in range(N_dir):
        directory_i = directories[i]
        directory = working_directory + "/" + directory_i
        files = os.listdir(directory)
        if os.path.exists(directory + '/sigma=15.000'):
            Z_r = np.loadtxt(directory + '/sigma=15.000/gamma=0.150/field_real.txt', delimiter=',')
            Z_i = np.loadtxt(directory + '/sigma=15.000/gamma=0.150/field_img.txt', delimiter=',')
            X = np.loadtxt(directory + '/sigma=15.000/gamma=0.150/X.txt', delimiter=',')
            T = np.loadtxt(directory + '/sigma=15.000/gamma=0.150/T.txt', delimiter=',')
            params = np.loadtxt(directory + '/sigma=15.000/gamma=0.150/parameters.txt', delimiter=',')
            sigma = 15#*sigmas[i]

            U_complex = Z_r + 1j * Z_i
            modulo = np.absolute(U_complex)
            arg = np.angle(U_complex)
            arg = (2 * np.pi + arg) * (arg < 0) + arg * (arg > 0)
            final_phase = arg[-1, :]
            final_amplitude = modulo[-1, :]
            forcing = params[2] * np.exp(- X ** 2 / (2 * sigma ** 2))
            delta_homo = - params[4] + np.sqrt(params[2] ** 2 + params[3] ** 2)


            def soliton(x, n, delta, x0):
                return n * np.sqrt(2 * delta) * (1 / (np.cosh(np.sqrt(delta) * (x - x0))))


            popt, pcov = curve_fit(soliton, X, final_amplitude, bounds=([0, 0, -30], [10, 10, 30]))
            parameter = (2 * np.arccosh(2) * np.sqrt(params[0] / delta_homo)) / (2 * np.sqrt(2 * np.log(2)) * sigma)
            print(parameter)
