from back_process import *


if __name__ == "__main__":
    disco = 'C'
    initial_dir_data = str(disco) + ':/Users/mnustes_science/PT_fluids/mnustes_science/simulation_data/simulation_data'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    y = np.loadtxt(directory + '/lyapunov_espectrum.txt', delimiter=',')
    plt.plot(y)
    plt.show()
    plt.close()
    N = len(y[:, 0])
    plt.plot(np.flip(np.sort(np.mean(y[int(0.2*N):], axis=0))))
    plt.xlim([0,200])
    plt.show()
    plt.close()

