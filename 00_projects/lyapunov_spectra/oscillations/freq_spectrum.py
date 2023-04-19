from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    disco = 'E:/'
    initial_dir_data = str(disco) + 'Users/mnustes_science/PT_fluids/mnustes_science/simulation_data/simulation_data'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')

    Z_r = np.loadtxt(directory + '/field_real_filt.txt', delimiter=',')
    Z_i = np.loadtxt(directory + '/field_img_filt.txt', delimiter=',')
    X = np.loadtxt(directory + '/X.txt', delimiter=',')
    T = np.loadtxt(directory + '/T.txt', delimiter=',')

    Nx = len(X)
    Nt = len(T)
    ti = int(0.1 * Nt)
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
    pcm = plt.pcolormesh(X, xf_i, np.transpose(Z_r_fft), cmap="jet", shading='auto')
    cbar = plt.colorbar(pcm)
    cbar.set_label('$A_R(x, \omega)$', rotation=0, size=25, labelpad=-27, y=1.11)

    # put the major ticks at the middle of each cell
    plt.xlim([X[0], X[-1]])
    plt.xlabel('$x$', size='25')
    plt.xticks(fontsize=20)

    plt.ylabel('$\omega$', size='25')
    plt.yticks(fontsize=20)
    plt.ylim([0, 0.2])

    plt.grid(linestyle='--', alpha=0.2, color='k')
    cbar.ax.tick_params(labelsize=15)

    # labels
    plt.tight_layout()
    plt.savefig(directory + '/field_real_FFT_filt.png', dpi=300)
    plt.close()