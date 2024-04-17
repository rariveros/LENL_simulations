from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    disco = 'E:/'
    initial_dir_data = str(disco) + "mnustes_science/experimental_data/soliton_control"
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    directories = [name for name in os.listdir(working_directory) if os.path.isdir(os.path.join(working_directory, name))]

    datafile = "C:/mnustes_science/simulation_data/FD/PDNLS_extended_PT/extras/dimensional/analysis"
    centers = []
    nus = []
    n = 0
    for directory in directories:
        if os.path.exists(working_directory + "/" + directory + '/centers.txt') and os.path.exists(working_directory + "/" + directory + '/nus.txt'):
            center_i = np.loadtxt(working_directory + "/" + directory + '/centers.txt', delimiter=',')
            nu_i = np.loadtxt(working_directory + "/" + directory + '/nus.txt', delimiter=',')
            centers.append(np.abs(center_i))
            nus.append(nu_i)
    nu = np.mean(nus, axis=0)
    center = np.mean(centers, axis=0)
    center_std = np.std(centers, axis=0)
    plt.errorbar(nu, center, yerr=center_std,  c="k", fmt='o')
    plt.show()