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
    data_out = "D:/mnustes_science/simulation_data/FD/PDNLS_extended_PT/oscillations_dnu/data_out"
    index_01 = []
    index_02 = []
    parameters = []

    for directory in directories:
        print("#############   " + directory + "   #############")
        Z_r = np.loadtxt(working_directory + '/' + directory + '/sigma=6.000/gamma=0.280/dist=20.000/field_real.txt', delimiter=',')
        Z_i = np.loadtxt(working_directory + '/' + directory + '/sigma=6.000/gamma=0.280/dist=20.000/field_img.txt', delimiter=',')
        X = np.loadtxt(working_directory + '/' + directory + '/sigma=6.000/gamma=0.280/dist=20.000/X.txt', delimiter=',')
        T = np.loadtxt(working_directory + '/' + directory + '/sigma=6.000/gamma=0.280/dist=20.000/T.txt', delimiter=',')
        params = np.loadtxt(working_directory + '/' + directory + '/sigma=6.000/gamma=0.280/dist=20.000/parameters.txt', delimiter=',')
        center_L = - 0.5 * params[3]
        center_R = 0.5 * params[3]

        Nx = len(X)
        Nt = len(T)

        I_L = np.where(X == center_L)[0][0]
        I_R = np.where(X == center_R)[0][0]

        ti, tf = int(0.5 * Nt), Nt
        Z_r = Z_r[ti:tf, :]
        Z_i = Z_i[ti:tf, :]
        T = T[ti:tf] - T[ti]
        dt = T[1] - T[0]

        Nx = len(X)
        Nt = len(T)

        Z = Z_r + 1j * Z_i
        Z_conj = Z_r - 1j * Z_i
        Z_modulo = np.absolute(Z)

        CCF = []
        CCF_flipped = []

        for j in range(Nx):
            CCF_flipped_j = np.correlate(np.flip(Z, axis=1)[:, j], Z_conj[:, j] ,  "full")
            CCF_j = np.correlate(Z[:, j], Z_conj[:, j], "full")
            CCF_flipped.append(CCF_flipped_j)
            CCF.append(CCF_j)
        CCF = np.array(np.abs(CCF))
        CCF_flipped = np.array(np.abs(CCF_flipped))
        tau = np.arange(-T[-1], T[-1], dt)

        middle_tau = int(len(tau) / 2)
        middle_x = int(len(X) / 2)
        S_L = np.transpose(CCF)[middle_tau, :middle_x]
        S_R = np.transpose(CCF)[middle_tau, middle_x:]
        X_L = X[:middle_x]
        X_R = X[middle_x:]

        S_L_integred = float(integrate.simpson(S_L, X_L))
        S_R_integred = float(integrate.simpson(S_R, X_R))
        delta_S = S_L - np.flip(S_R)

        int_parameter = np.abs((S_L_integred - S_R_integred) / (S_L_integred + S_R_integred))   #va de cero a uno
        desvesta_parameter = np.mean(delta_S ** 2) #no normalizado
        index_01.append(int_parameter)
        index_02.append(desvesta_parameter)
        parameters.append(params[5])

    if not os.path.exists(data_out):
        os.makedirs(data_out)

    plt.scatter(parameters, index_01)
    plt.xlabel('$\\nu$', size='20')
    plt.ylabel('$\\xi^{1}_s$', size='20')
    #plt.ylim(0.05, 0.125)
    plt.grid(linestyle='--', alpha=0.5)
    plt.savefig(data_out + '//index_01.png', dpi=300)
    plt.close()

    plt.scatter(parameters, index_02)
    plt.xlabel('$\\nu$', size='20')
    plt.ylabel('$\\xi^{2}_s$', size='20')
    #plt.ylim(0.05, 0.125)
    plt.grid(linestyle='--', alpha=0.5)
    plt.savefig(data_out + '//index_02.png', dpi=300)
    plt.close()

