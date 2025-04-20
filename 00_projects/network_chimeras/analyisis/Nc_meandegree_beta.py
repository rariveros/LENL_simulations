from back_process import *

if __name__ == '__main__':
    y = [0, 0, 0.05, 0.16, 0.30, 0.49, 0.62, 0.68, 0.78, 0.83, 0.92, 0.96, 0.93, 0.94, 0.97, 0.95]
    x = [15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]

    fig, ax = plt.subplots(1, 1, figsize=(3, 3))
    ax.scatter(x, y, color="k")
    ax.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
    ax.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=True)
    ax.set_ylabel("$\\frac{N_c}{N}$", fontsize=20)
    ax.set_xlabel("$\\langle k \\rangle$", fontsize=20)
    ax.text(21, 0.1, "$\kappa = 1.5 \\times 10^{-2}$", fontsize=15)
    ax.text(21, 0., "$N = 501$", fontsize=15)
    plt.tight_layout()
    plt.savefig('Fig03.png', dpi=200)

    # HIPOTESIS 0: en osciladores de duffing, para tales parametros (referencias), coexisten quimeras (caos) y estados sincronizados.
    # HIPOTESIS 1: k efectivo atrae a la quimera, incluso en redes de baja clusterización hay correlación entre centralidad y caos.
    # HIPOTESIS 2: el estado estacionario final es independiente (estadisticamente) del pinch de quimera inicial

    # CARACTERIZACIÓN 0: espectro de lyapunov, correlaciones con fase/modulo, correlaciones con centralidad e implicancias de esto.
    # CARACTERIZACIÓN 1: variar k y mean degree para erdos renyi, definir ancho efectivo (delta) en mean degree y graficar versus k (busqueda de parametro de control).
    # (es una transición? que tipo?)
    # CARACTERIZACIÓN 2: probar que para una red, el estado estacionario es estadisticamente independiente del pinch (como visualizar?)
    # APLICACIÓN 1: red de alta clusterización, identificar clusters con alto delta, proponer para clasificación de comunidades (acá el trabajo se puede diseccionar).