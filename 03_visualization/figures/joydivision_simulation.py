from directories import *
from back_process import *

from matplotlib.animation import FuncAnimation, PillowWriter

if __name__ == "__main__":
    disco = 'C'
    initial_dir_data = str(disco) + ':/Users/mnustes_science/PT_fluids/mnustes_scienc'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    A = np.loadtxt(directory + '/A_real.txt', delimiter=',', dtype=complex)
    C_R = np.loadtxt(directory + '/field_real.txt', delimiter=',', dtype=complex)
    C_I = np.loadtxt(directory + '/field_img.txt', delimiter=',', dtype=complex)
    C_mod = np.sqrt(C_R ** 2 + C_I ** 2)
    X = np.loadtxt(directory + '/X.txt', delimiter=',')
    T = np.loadtxt(directory + '/T.txt', delimiter=',')

    Nt = len(T)
    print(Nt)
    t_200 = int(Nt / 5)
    dt = int(t_200 / 200)
    t_init = t_200 + 42 * dt
    N = 600

    plt.plot(X, np.real(A[int(t_init), :])/2, color="r")
    plt.plot(X, np.real(C_mod[int(t_init), :]), color="k", linestyle="--")
    plt.plot(X, -np.real(C_mod[int(t_init), :]), color="k", linestyle="--")
    px = 1 / plt.rcParams['figure.dpi']
    #fig1, ax1 = plt.subplots(figsize=(1200 * px, 1400 * px))
    #for n in range(0, Nt, 1):
    #    ax1.plot(X, np.real(C_mod[int(t_init + n * dt), :]) + n * 2, color=(1 - (n / Nt), 0, n / Nt, 0.5))

    fig2, ax2 = plt.subplots(figsize=(1200 * px, 1400 * px))
    for n in range(0, N, 13):
        ax2.plot(X, np.abs(A[int(t_init + n * dt), :])/2 + n * 0.01, color=(1 - (n / N), 0, n / N, 0.95), lw=3, zorder=n)
        #ax2.plot(X, np.real(C_mod[int(t_init + n * dt), :]) + n * 0.5, color="k", linestyle="--", lw=2, alpha=0.85, zorder=n)
        ax2.set_xlim([-35, 35])
        plt.yticks([])
        plt.xticks([-30, -20, -10, 0, 10, 20, 30])
        ax2.set_xlabel('$x$', fontsize=50)
        ax2.tick_params(labelsize=45)
        #plt.savefig('E:/mnustes_science/experimental_data/faraday_drift_03/dvelocities_info/velocity_processed/fit', dpi=300)

    #fig3, ax3 = plt.subplots(figsize=(1200 * px, 1400 * px))
    #for n in range(0, Nt, 1):
    #    ax3.plot(X, np.abs(np.real(A[int(t_init + n * dt), :])) + n * 2, color=(1 - (n / Nt), 0, n / Nt, 0.5), zorder=n)
    #plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(directory + '/C_module.png', dpi=300)
    plt.close()

