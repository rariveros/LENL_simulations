import matplotlib.pyplot as plt
from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    frequencies = []
    disc = "aD:"
    initial_dir_data = str(disc) + 'mnustes_science/simulation_data/FD'
    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=initial_dir_data, title='Elección de carpeta')
    VR = np.loadtxt(directory + '/v_R.txt', delimiter=',')
    VL = np.loadtxt(directory + '/v_L.txt', delimiter=',')
    AR = np.loadtxt(directory + '/A_R.txt', delimiter=',')
    AL = np.loadtxt(directory + '/A_L.txt', delimiter=',')
    TR = np.loadtxt(directory + '/T_R.txt', delimiter=',')
    TL = np.loadtxt(directory + '/T_L.txt', delimiter=',')
    params = np.loadtxt(directory + '/parameter.txt', delimiter=',')


    V = np.array([np.abs(VR), np.abs(VL)])
    T = np.array([np.abs(TR), np.abs(TL)])
    A = np.array([np.array(AR), np.array(AL)])

    # Ordenar filas para que V[0, :] sea siempre el menor valor entre VR y VL
    i_min = np.argmin(V, axis=0)
    i_max = 1 - i_min

    cols = np.arange(V.shape[1])

    V_sorted = np.vstack([V[i_min, cols], V[i_max, cols]])
    T_sorted = np.vstack([T[i_min, cols], T[i_max, cols]])
    A_sorted = np.vstack([A[i_min, cols], A[i_max, cols]])

    # Ordenar por parámetro
    j_order = np.argsort(params)
    params_sorted = params[j_order]
    V_sorted = V_sorted[:, j_order]
    T_sorted = T_sorted[:, j_order]
    A_sorted = A_sorted[:, j_order]

    beta = 0.004811649356064012

    # Rango de parámetros
    param_min = 0.12
    param_max = 0.23

    # Tamaños para presentación
    label_fontsize = 28
    tick_fontsize = 20
    title_fontsize = 18
    point_size = 60

    fig, (ax01, ax02, ax04, ax03) = plt.subplots(1, 4, figsize=(21, 4))

    # --- Gráfico superior izquierdo ---
    ax01.scatter(params_sorted, V_sorted[0], label='$\\textrm{Right}$', color='green', s=point_size)
    ax01.scatter(params_sorted, V_sorted[1], label='$\\textrm{Left}$', color='purple', s=point_size)
    ax01.set_ylabel('$\langle v \\rangle$', fontsize=label_fontsize)
    ax01.set_xlim(param_min, param_max)
    #ax01.set_ylim(-0.025, 0.25)
    ax01.tick_params(labelsize=tick_fontsize, direction='in')
    ax01.grid(True, alpha=0.3)
    ax01.hlines(0, param_min, param_max, colors='black', linestyles='-')

    # --- Gráfico superior derecho ---
    ax02.scatter(params_sorted, A_sorted[0] / np.sqrt(beta), color='green', s=point_size)
    ax02.scatter(params_sorted, A_sorted[1] / np.sqrt(beta), color='purple', s=point_size)
    ax02.set_xlim(param_min, param_max)
    ax02.tick_params(labelsize=tick_fontsize, direction='in')
    #ax02.set_ylabel('$|\psi|_{\\textrm{max}}$', fontsize=label_fontsize, labelpad=20)
    ax02.grid(True, alpha=0.3)
    ax02.hlines(0, param_min, param_max, colors='black', linestyles='-')

    # --- Gráfico inferior derecho ---
    ax03.scatter(params_sorted, V_sorted[0] * T_sorted[0], color='green', s=point_size)
    ax03.scatter(params_sorted, V_sorted[1] * T_sorted[1], color='purple', s=point_size)
    ax03.set_xlabel('$\gamma_0$', fontsize=label_fontsize)
    ax03.set_ylabel('$\langle \\lambda \\rangle$', fontsize=label_fontsize, labelpad=20)
    ax03.set_xlim(param_min, param_max)
    ax03.tick_params(labelsize=tick_fontsize, direction='in')
    ax03.grid(True, alpha=0.3)
    ax03.hlines(0, param_min, param_max, colors='black', linestyles='-')

    # --- Gráfico inferior izquierdo ---
    ax04.scatter(params_sorted, T_sorted[0], color='green',s=point_size)
    ax04.scatter(params_sorted, T_sorted[1], color='purple', s=point_size)
    ax04.set_xlabel('$\gamma_0$', fontsize=label_fontsize)
    ax04.set_ylabel('$\langle T \\rangle$', fontsize=label_fontsize)
    ax04.set_xlim(param_min, param_max)
    ax04.tick_params(labelsize=tick_fontsize, direction='in')
    ax04.grid(True, alpha=0.3)
    ax04.hlines(0, param_min, param_max, colors='black', linestyles='-')

    # Layout ajustado
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.1, wspace=0.25)
    plt.savefig(directory + "/characterization.png", dpi=300)