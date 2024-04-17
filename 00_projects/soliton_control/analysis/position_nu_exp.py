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
    if not os.path.exists(datafile):
        os.makedirs(datafile)
    n = 0
    for directory in directories:
        print("#############   " + directory + "   #############")
        Z = np.loadtxt(working_directory + "/" + directory + '/Z_mm.txt', delimiter=',')
        T = np.loadtxt(working_directory + "/" + directory + '/T_s.txt', delimiter=',')
        X = np.loadtxt(working_directory + "/" + directory + '/X_mm.txt', delimiter=',')
        marks = np.loadtxt(working_directory + "/" + directory + '/IL_mm.txt', delimiter=',')
        print(marks)
        x_0 = np.mean(marks)
        print(x_0)
        X = X - x_0
        file_name = os.path.basename(directory)
        name_list = file_name.split("_")
        notation = 'af'
        if notation == 'fa':
            f_i = float(name_list[0].split("=")[-1])
            a = float(name_list[1].split("=")[-1])
        elif notation == 'af':
            a = float(name_list[0].split("=")[-1])
            f_i = float(name_list[1].split("=")[-1])
        elif notation == "respuesta":
            a = 0
            f_i = 7
        d = 20
        alpha, beta, nu, gamma, f_0 = fluid_pdnls_parameters(f_i, a, d)
        #[alpha, beta, mu_0, nu, sigma, gamma_0]
        beta = 1
        Z_mean = np.mean(np.abs(Z), axis=0)

        Nx = len(X)
        Nt = len(T)

        dt = T[1] - T[0]

        Nx = len(X)
        Nt = len(T)

        Z_last = Z_mean

        J_max = np.argmax(Z_last)
        delta_wind = 50
        window_L = J_max - delta_wind
        window_R = J_max + delta_wind

        def quadratic(x, a, b, c):
            return a * x ** 2 + b * x + c

        popt, pcov = curve_fit(quadratic, X[window_L:window_R], Z_last[window_L:window_R])
        #plt.scatter(X[window_L:window_R], Z_last[window_L:window_R])
        #plt.plot(X[window_L:window_R], quadratic(X[window_L:window_R], *popt))
        #plt.show()
        #plt.close()
        #plt.plot(X, Z_last + n)
        [a, b, c] = popt
        center = - b / (2 * a)
        centers.append(center)
        nus.append(nu)
        n = n + 0.5
    #plt.show()
    #plt.close()
    np.savetxt(working_directory + '/centers.txt', centers, delimiter=',')
    np.savetxt(working_directory + '/nus.txt', nus, delimiter=',')
    #plt.scatter(nus, np.abs(centers), c="k", zorder=10)
    #plt.xlabel("$\\nu$", fontsize=20)
    #plt.ylabel("$\Delta x$", fontsize=20)
    #plt.grid(alpha=0.2, zorder=0)
    #plt.show()

