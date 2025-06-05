from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    frequencies = []
    disc = "D:"
    initial_dir_data = str(disc) + 'mnustes_science/simulation_data/FD'
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    directories = [name for name in os.listdir(working_directory) if os.path.isdir(os.path.join(working_directory, name))]
    for directory_01 in directories:
        dir_01 = working_directory + "/" + directory_01
        directories_01 = [name for name in os.listdir(dir_01) if os.path.isdir(os.path.join(dir_01, name))]
        frequencies = []
        omegas = []
        ks = []
        for directory_02 in directories_01:
            dir_02 = working_directory + "/" + directory_01 + "/" + directory_02
            params = np.loadtxt(dir_02 + '/parameters.txt', delimiter=',')
            X = np.loadtxt(dir_02 + '/X.txt', delimiter=',')
            T = np.loadtxt(dir_02 + '/T.txt', delimiter=',')
            ZR = np.loadtxt(dir_02 + '/field_real.txt', delimiter=',')
            ZI = np.loadtxt(dir_02 + '/field_img.txt', delimiter=',')

            Nt = len(T)
            I_inicial = int(0.6 * Nt)
            I_final = -1
            dI = 10
            T = T[I_inicial:I_final:dI]
            ZR = ZR[I_inicial:I_final:dI, :]
            ZI = ZI[I_inicial:I_final:dI, :]

            save_directory = "D:/mnustes_science/simulation_data/FD/bi_localized_drift/datos_fernando/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.1000/sigma=16.000" + "/" + directory_01
            if not os.path.exists(save_directory):
                os.makedirs(save_directory)
            np.savetxt(save_directory + '/field_real.txt', ZR, delimiter=',')
            np.savetxt(save_directory + '/field_img.txt', ZI, delimiter=',')
            np.savetxt(save_directory + '/parameters.txt', params, delimiter=',')
            np.savetxt(save_directory + '/X.txt', X, delimiter=',')
            np.savetxt(save_directory + '/T.txt', T, delimiter=',')