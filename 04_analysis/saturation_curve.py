from directories import *
from back_process import *

if __name__ == "__main__":
    disco = 'E'
    initial_dir_data = str(disco) + ':/mnustes_science/simulation_data'
    save_directory = initial_dir_data + '/tesis_rafael/pattern_saturation'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    directories = [name for name in os.listdir(directory) if os.path.isdir(os.path.join(directory, name))]
    for i in range(len(directories)):
        current = directory + '/' + directories[i]
        currents = [name for name in os.listdir(current) if os.path.isdir(os.path.join(current, name))]
        current = current + '/' + currents[-1]
        currents = [name for name in os.listdir(current) if os.path.isdir(os.path.join(current, name))]
        current = current + '/' + currents[-1]
        print(current)
        nu_str = current.split('/')[-2].split('=')[-1]
        gamma_str = current.split('/')[-3].split('=')[-1]
        mu_str = current.split('/')[-4].split('=')[-1]

        print('mu = ' + mu_str + '    gamma = ' + gamma_str + '    nu = ' + nu_str)

        nu = float(nu_str)
        gamma = float(gamma_str)
        mu = float(mu_str)

        sigma = 10


        X = np.loadtxt(current + '/X.txt', delimiter=',')
        T = np.loadtxt(current + '/T.txt', delimiter=',')
        A_r = np.loadtxt(current + '/field_real.txt', delimiter=',')
        A_i = np.loadtxt(current + '/field_img.txt', delimiter=',')

        Nx = len(X)
        Nt = len(T)
        B = np.zeros((Nt, Nx))
        A = A_r + 1j * A_i
        module = np.absolute(A)
        phase = np.angle(A)
        phase = (2 * np.pi + phase) * (phase < 0) + phase * (phase > 0)

        signal_r = hilbert(A_r)
        B_r = np.abs(signal_r)
        B_r = np.mean(B_r, axis=1)

        signal_i = hilbert(A_i)
        B_i = np.abs(signal_i)
        B_i = np.amax(B_i, axis=1)

        B = B_r + 1j * B_i
        B_0 = (2 * mu / 9) ** 0.25 * (gamma - (mu + (1 / sigma) * np.sqrt(nu))) ** 0.25

        if not os.path.exists(save_directory):
            os.makedirs(save_directory)

        np.savetxt(save_directory + '/B_' + 'mu=' + mu_str + '_gamma=' + gamma_str + '_nu=' + nu_str + '.txt', B, delimiter=',')
        plt.plot(T, np.abs(B))
        plt.hlines(B_0, T[0], T[-1])
        plt.savefig(save_directory + '/B_' + 'mu=' + mu_str + '_gamma=' + gamma_str + '_nu=' + nu_str + '.png', dpi=300)
        plt.close()