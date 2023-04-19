import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *
from directories_lyap import *

if __name__ == '__main__':
    disco = 'E:/'
    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FD' + main_directory + subdir_exponent
    root = tk.Tk()
    root.withdraw()
    working_directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    D_ky = np.loadtxt(working_directory + '/KY_dim_sigma.txt', delimiter=',')
    sigma = np.loadtxt(working_directory + '/sigmas.txt', delimiter=',')

    ### Ajuste de curva ###
    def linear_fit(x, m, n):
        return m * x + n

    def sqr_root(x, A, c):
        x = x.astype('complex128')
        return A * np.abs((x - c) ** 0.5)


    popt, pcov = curve_fit(linear_fit, sigma, D_ky)

    m = popt[0]
    n = popt[1]

    print('m='+str(m))
    print('n=' + str(n))

    sigma_fit = np.arange(0, 22, 0.1)
    D_ky_fit = linear_fit(sigma_fit, *popt)

    residuals = D_ky - linear_fit(sigma, *popt)
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((D_ky - np.mean(D_ky)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)
    print(r_squared)


    px = 1/plt.rcParams['figure.dpi']
    fig, ax1 = plt.subplots(figsize=(900*px, 800*px))
    textstr = '\n'.join((r'$a=%.3f$' % (m,), r'$b=%.3f$' % (n,), r'$R^2=%.3f$' % (r_squared,)))
    textstr_02 = "$D_{KY} = a \sigma_i + b$"
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax1.text(0.05, 0.95, textstr, transform=ax1.transAxes, fontsize=20, verticalalignment='top', bbox=props)
    ax1.set_xlabel('$\sigma_i$', fontsize=30)
    ax1.set_ylabel('$D_{KY}$', fontsize=30)
    ax1.tick_params(labelsize=20)
    ax1.plot(sigma_fit, D_ky_fit, '--', linewidth='3', alpha=0.8, c='r', zorder=0)
    ax1.scatter(sigma, D_ky, color="k", zorder=1)
    ax1.set_ylim([0, 50])
    ax1.set_xlim([5, 21])
    ax1.grid(alpha=0.2, color="k", zorder=0)
    plt.tight_layout()
    plt.savefig(working_directory + "/D_ky_fit.png", dpi=300)