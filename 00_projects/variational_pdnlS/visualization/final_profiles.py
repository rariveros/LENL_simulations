from back_process import *


if __name__ == '__main__':
    #disco = 'C:/'
    #initial_dir_data = str(disco) + ':/mnustes_science/simulation_data'
    #root = tk.Tk()
    #root.withdraw()
    #directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    directory = "C:/mnustes_science/simulation_data/FD/PT_extended/variando_distancia/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.0180/sigma=6.000/dist=32.000"
    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    X = np.loadtxt(directory + '/X.txt', delimiter=',')
    T = np.loadtxt(directory + '/T.txt', delimiter=',')

    Z = Z_r + 1j * Z_i
    Z_conj = Z_r - 1j * Z_i
    Z_modulo = np.absolute(Z)
    plt.plot(X, Z_r[-1, :],color="b")
    plt.plot(X, Z_i[-1, :],color="r")
    plt.plot(X, np.abs(Z[-1, :]), color="k")
    plt.hlines(0., -120, 120, colors="k", alpha=0.4)
    plt.grid(alpha=0.2)
    plt.show()