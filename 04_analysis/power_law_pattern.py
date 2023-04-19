from directories import *
from back_process import *

if __name__ == "__main__":
    disco = 'E'
    initial_dir_data = str(disco) + ':/mnustes_science/simulation_data'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    dir_list = os.listdir(directory)
    Bs = []
    deltas = []
    for file in dir_list:
        parameters = file.split('_')
        nu_str = parameters[-1].split('=')[-1][0:-4]
        gamma_str = parameters[-2].split('=')[-1]
        mu_str = parameters[-3].split('=')[-1]

        nu = float(nu_str)
        gamma = float(gamma_str)
        mu = float(mu_str)

        B = np.loadtxt(directory + '/' + file, delimiter=',', dtype=np.complex_)
        B_module = np.abs(B)
        Bs.append(B_module[-1])
        delta = gamma - mu
        deltas.append(delta)
    np.savetxt('E:/mnustes_science/simulation_data/tesis_rafael/pattern_powerlaw/Bs.txt', Bs, delimiter=',')
    np.savetxt('E:/mnustes_science/simulation_data/tesis_rafael/pattern_powerlaw/deltas.txt', deltas, delimiter=',')
    plt.scatter(deltas, Bs, c='k', zorder=3)
    plt.rc('axes', axisbelow=True)
    plt.grid(zorder=0)
    plt.xlabel('$\delta$', rotation=0, size=20)
    plt.ylabel('$B_0$', rotation=0, size=20, labelpad=15)
    plt.savefig('power_law.png', dpi=300)
    plt.close()


