import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    disco = 'D:/'
    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FDa'
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    directories = [name for name in os.listdir(working_directory) if os.path.isdir(os.path.join(working_directory, name))]

    datafile = "C:/mnustes_science/simulation_data/FD/PDNLS_extended_PT/extras/dimensional/analysis"
    centers = []
    nus = []
    gammas = ["0.180"]
    sigmas = ["15.000", "17.500", "20.000", "22.500", "25.000", "27.500", "30.000"]
    if not os.path.exists(datafile):
        os.makedirs(datafile)
    for gamma in gammas:
        for i in range(len(sigmas)):
            sigma = sigmas[i]
            deltas = []
            amps = []
            anchos = []
            for directory in directories:
                if os.path.exists(working_directory + "/" + directory + '/sigma=' + sigma + '/gamma=' + gamma):
                    data = np.loadtxt(
                        working_directory + "/" + directory + '/sigma=' + sigma + '/gamma=' + gamma + '/data.txt',
                        delimiter=',')
                    [alpha, beta, gamma_0, mu_0, nu, sig] = np.loadtxt(
                        working_directory + "/" + directory + '/sigma=' + sigma + '/gamma=' + gamma + '/parameters.txt',
                        delimiter=',')
                    # [alpha, beta, gamma_0, mu_0, nu, sigma]
                    # [center, delta, ancho, amplitud]
                    center = data[0]
                    deltas.append(np.abs(center))
                    plt.scatter(nu, np.abs(center), c="k")
    plt.show()
