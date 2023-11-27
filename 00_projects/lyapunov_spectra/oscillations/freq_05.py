from functions import *
from back_process import *
from time_integrators import *
import scipy
from scipy import interpolate

if __name__ == '__main__':
    disco = 'E:/'
    initial_dir_data = str(disco) + 'mnustes_science/simulation_data/FDa'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    freq = np.loadtxt(directory + '/freq.txt', delimiter=',')
    param = np.loadtxt(directory + '/param.txt', delimiter=',')
    result = scipy.stats.linregress(param[2:], freq[2:])
    print(result.slope, result.intercept)
    x = np.arange(0.158, 0.22, 0.001)
    y = result.slope * x + result.intercept
    plt.plot(x, y, c="r", linestyle="--", linewidth=3)
    plt.scatter(param, freq, c="k", zorder=10, s=50)
    plt.grid(alpha=0.5, zorder=1)
    plt.ylim([0, 0.35])
    plt.xlim([0.148, 0.2])
    plt.xlabel("$\gamma_0$", fontsize=25)
    plt.ylabel("$f\ \\textrm{(Hz)}$", fontsize=25)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tight_layout()
    plt.savefig('freq_example.png', dpi=300)
    plt.close()