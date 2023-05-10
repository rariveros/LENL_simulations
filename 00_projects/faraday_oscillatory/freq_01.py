from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    disco = 'D:/'
    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FD/PDNLS_oscillatory/alpha=5.721/beta=1.000/mu=0.100/nu=0.056/sigma=60.000'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
    X = np.loadtxt(directory + '/X.txt', delimiter=',')
    T = np.loadtxt(directory + '/T.txt', delimiter=',')

    Nx = len(X)
    Nt = len(T)
    ti = int(0.75 * Nt)
    tf = -1

    Z_r = Z_r[ti:tf, :]
    T = T[ti:tf]
    Nt = len(T)
    dt = T[1] - T[0]
    Z_r_fft = []
    for i in range(Nx):
        yf_i = fft(Z_r[:, i])
        yf_i = 2.0 / Nt * np.abs(yf_i[0:Nt // 2])
        xf_i = fftfreq(Nt, dt)[:Nt // 2]
        Z_r_fft.append(yf_i.tolist())

    Z_r_fft = np.array(Z_r_fft)
    fig, ax = plt.subplots()

    # legend
    pcm = plt.pcolormesh(X, xf_i[0:], np.log(np.transpose(Z_r_fft[:, 0:])), vmin=-10, vmax=0, cmap="jet", shading='auto')
    cbar = plt.colorbar(pcm)
    cbar.set_label('$A_R(x, \omega)$', rotation=0, size=25, labelpad=-27, y=1.11)

    # put the major ticks at the middle of each cell
    plt.xlim([X[0], X[-1]])
    plt.xlabel('$x$', size='25')
    plt.xticks(fontsize=20)

    plt.ylabel('$\omega$', size='25')
    plt.yticks(fontsize=20)
    plt.ylim([0, 1])

    plt.grid(linestyle='--', alpha=0.2, color='k')
    cbar.ax.tick_params(labelsize=15)

    # labels
    plt.tight_layout()
    plt.savefig(directory + '/field_real_FFT.png', dpi=300)
    plt.close()