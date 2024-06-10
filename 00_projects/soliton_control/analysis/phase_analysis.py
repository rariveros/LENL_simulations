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
    colors = (np.array([10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30]) - 10) / 20
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
                    print("#############   " + directory + "   #############")
                    Z_r = np.loadtxt(working_directory + "/" + directory + '/sigma=' + sigma + '/gamma=' + gamma + '/field_real.txt', delimiter=',')
                    Z_i = np.loadtxt(working_directory + "/" + directory + '/sigma=' + sigma + '/gamma=' + gamma + '/field_img.txt', delimiter=',')
                    X = np.loadtxt(working_directory + "/" + directory + '/sigma=' + sigma + '/gamma=' + gamma + '/X.txt', delimiter=',')
                    [alpha, beta, gamma_0, mu_0, nu, sig] = np.loadtxt(working_directory + "/" + directory + '/sigma=' + sigma + '/gamma=' + gamma + '/parameters.txt', delimiter=',')
                    Z = Z_r[-1, :] + 1j * Z_i[-1, :]
                    arg = np.angle(Z)
                    gam = gamma_0 * np.exp(- X ** 2 / (2 * sig ** 2))
                    c = gam * np.cos(2 * arg)
                    s = gam * np.sin(2 * arg)
                    #[alpha, beta, gamma_0, mu_0, nu, sigma]
                    #[center, delta, ancho, amplitud]
                    plt.plot(X, c, color="b")
                    plt.plot(X, s, color="r")
                    plt.show()
                    plt.close()
