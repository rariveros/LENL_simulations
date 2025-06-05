from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    frequencies = []
    disc = "aD:"
    initial_dir_data = str(disc) + 'mnustes_science/simulation_data/FD'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    distances = np.loadtxt(directory + '/analysis/dists.txt', delimiter=',')
    T = np.loadtxt(directory +  '/analysis/t_grid.txt', delimiter=',')
    UR_Rs = np.loadtxt(directory + '/analysis/UR_Rs.txt', delimiter=',')
    UR_Is = np.loadtxt(directory + '/analysis/UR_Is.txt', delimiter=',')
    UL_Rs = np.loadtxt(directory + '/analysis/UL_Rs.txt', delimiter=',')
    UL_Is = np.loadtxt(directory + '/analysis/UL_Is.txt', delimiter=',')
    save_directory = directory + '/analysis'

    for i in range(len(distances)):
        plt.plot(T, UR_Rs[i])
        plt.show()
