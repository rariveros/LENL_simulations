from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':

    def function(phi1, phi2, k, params):
        [mu, gamma, nu] = params
        R = np.sqrt(((gamma + mu) / gamma) * (-nu + np.sqrt(gamma ** 2 - mu ** 2)))
        P = np.sqrt(R ** 2 * (gamma - mu) / (gamma + mu))

        NL01a = P ** 2 * (np.cos(phi1) * np.cos(phi2) ** 3 + np.sin(phi1) * np.sin(phi2) ** 3)
        NL02a = R ** 2 * (np.cos(phi1) * np.sin(phi1) ** 2 * np.cos(phi2) + np.cos(phi1) ** 2 * np.sin(phi1) * np.sin(phi2))
        NL01b = R ** 2 * (np.cos(phi2) * np.cos(phi1) ** 3 + np.sin(phi2) * np.sin(phi1) ** 3)
        NL02b = P ** 2 * (np.cos(phi2) * np.sin(phi2) ** 2 * np.cos(phi1) + np.cos(phi2) ** 2 * np.sin(phi2) * np.sin(phi1))

        F1 = - k - nu * (P / R) * np.cos(phi1 - phi2) - (P / R) * (NL01a + NL02a)
        F2 = - k - nu * (R / P) * np.cos(phi1 - phi2) - (R / P) * (NL01b + NL02b)
        return np.array([F1, F2])

    def jacobian(phi1, phi2, k, params):
        [mu, gamma, nu] = params
        R = np.sqrt(((gamma + mu) / gamma) * (-nu + np.sqrt(gamma ** 2 - mu ** 2)))
        P = np.sqrt(R ** 2 * (gamma - mu) / (gamma + mu))
        J11 = (-nu * (P / R) * np.sin(phi2 - phi1) + (P ** 3 / R) * (np.sin(phi1) * np.cos(phi2) ** 3 - np.cos(phi1) * np.sin(phi2) ** 3)
               - R * P * ((2 * np.sin(phi1) * np.cos(phi1) ** 2 - np.sin(phi1)) * np.cos(phi2) + (np.cos(phi1) ** 3 - 2 * np.sin(phi1) ** 2 * np.cos(phi1)) * np.sin(phi2) ** 3))
        J12 = (nu * (P / R) * np.sin(phi2 - phi1) + (3 * P ** 3 / R) * (np.cos(phi1) * np.cos(phi2) ** 2 * np.sin(phi2) - np.sin(phi1) * np.cos(phi2) * np.sin(phi2) ** 2)
               + R * P * (np.cos(phi1) * np.sin(phi1) ** 2 * np.sin(phi2) - np.cos(phi1) ** 2 * np.sin(phi1) * np.cos(phi2)))
        J21 = (-nu * (R / P) * np.sin(phi2 - phi1) + (3 * R ** 3 / P) * (np.cos(phi2) * np.cos(phi1) ** 2 * np.sin(phi1) - np.sin(phi2) * np.cos(phi1) * np.sin(phi1) ** 2)
               + R * P * (np.cos(phi2) * np.sin(phi2) ** 2 * np.sin(phi1) - np.cos(phi2) ** 2 * np.sin(phi2) * np.cos(phi1)))
        J22 = (nu * (R / P) * np.sin(phi2 - phi1) + (R ** 3 / P) * (np.sin(phi2) * np.cos(phi1) ** 3 - np.cos(phi2) * np.sin(phi1) ** 3)
               - R * P * ((2 * np.sin(phi2) * np.cos(phi2) ** 2 - np.sin(phi2)) * np.cos(phi1) + (np.cos(phi2) ** 3 - 2 * np.sin(phi2) ** 2 * np.cos(phi2)) * np.sin(phi1) ** 3))
        return np.array([[J11, J12], [J21, J22]])

    def H_fun(phi1, phi2, k, params, U_pred, dU_s):
        x_pred = U_pred[0]
        y_pred = U_pred[1]
        mu_pred = U_pred[2]

        dx_s = dU_s[0]
        dy_s = dU_s[1]
        dm_s = dU_s[2]

        [mu, gamma, nu] = params
        R = np.sqrt(((gamma + mu) / gamma) * (-nu + np.sqrt(gamma ** 2 - mu ** 2)))
        P = np.sqrt(R ** 2 * (gamma - mu) / (gamma + mu))

        NL01a = P ** 2 * (np.cos(phi1) * np.cos(phi2) ** 3 + np.sin(phi1) * np.sin(phi2) ** 3)
        NL02a = R ** 2 * (np.cos(phi1) * np.sin(phi1) ** 2 * np.cos(phi2) + np.cos(phi1) ** 2 * np.sin(phi1) * np.sin(phi2))
        NL01b = R ** 2 * (np.cos(phi2) * np.cos(phi1) ** 3 + np.sin(phi2) * np.sin(phi1) ** 3)
        NL02b = P ** 2 * (np.cos(phi2) * np.sin(phi2) ** 2 * np.cos(phi1) + np.cos(phi2) ** 2 * np.sin(phi2) * np.sin(phi1))

        F1 = - k - nu * (P / R) * np.cos(phi1 - phi2) - (P / R) * (NL01a + NL02a)
        F2 = - k - nu * (R / P) * np.cos(phi1 - phi2) - (R / P) * (NL01b + NL02b)
        H = (phi1 - x_pred) * dx_s + (phi2 - y_pred) * dy_s + (k - mu_pred) * dm_s
        return np.array([F1, F2, H])

    def monster_jacobian(phi1, phi2, k, params, dU_s):
        dx_s = dU_s[0]
        dy_s = dU_s[1]
        dm_s = dU_s[2]
        [mu, gamma, nu] = params

        R = np.sqrt(((gamma + mu) / gamma) * (-nu + np.sqrt(gamma ** 2 - mu ** 2)))
        P = np.sqrt(R ** 2 * (gamma - mu) / (gamma + mu))
        J11 = (-nu * (P / R) * np.sin(phi2 - phi1) + (P ** 3 / R) * (
                    np.sin(phi1) * np.cos(phi2) ** 3 - np.cos(phi1) * np.sin(phi2) ** 3)
               - R * P * ((2 * np.sin(phi1) * np.cos(phi1) ** 2 - np.sin(phi1)) * np.cos(phi2) + (
                            np.cos(phi1) ** 3 - 2 * np.sin(phi1) ** 2 * np.cos(phi1)) * np.sin(phi2) ** 3))
        J12 = (nu * (P / R) * np.sin(phi2 - phi1) + (3 * P ** 3 / R) * (
                    np.cos(phi1) * np.cos(phi2) ** 2 * np.sin(phi2) - np.sin(phi1) * np.cos(phi2) * np.sin(phi2) ** 2)
               + R * P * (np.cos(phi1) * np.sin(phi1) ** 2 * np.sin(phi2) - np.cos(phi1) ** 2 * np.sin(phi1) * np.cos(
                    phi2)))
        J13 = -1.0
        J21 = (-nu * (R / P) * np.sin(phi2 - phi1) + (3 * R ** 3 / P) * (
                    np.cos(phi2) * np.cos(phi1) ** 2 * np.sin(phi1) - np.sin(phi2) * np.cos(phi1) * np.sin(phi1) ** 2)
               + R * P * (np.cos(phi2) * np.sin(phi2) ** 2 * np.sin(phi1) - np.cos(phi2) ** 2 * np.sin(phi2) * np.cos(
                    phi1)))
        J22 = (nu * (R / P) * np.sin(phi2 - phi1) + (R ** 3 / P) * (
                    np.sin(phi2) * np.cos(phi1) ** 3 - np.cos(phi2) * np.sin(phi1) ** 3)
               - R * P * ((2 * np.sin(phi2) * np.cos(phi2) ** 2 - np.sin(phi2)) * np.cos(phi1) + (
                            np.cos(phi2) ** 3 - 2 * np.sin(phi2) ** 2 * np.cos(phi2)) * np.sin(phi1) ** 3))

        J23 = -1.0
        J31 = dx_s
        J32 = dy_s
        J33 = dm_s
        return np.array([[J11, J12, J13], [J21, J22, J23], [J31, J32, J33]])

    # --- semilla ---
    k = 0.03
    x0 = -0.4 * np.pi - 2 * np.pi
    y0 = 0.3 * np.pi - 2 * np.pi
    U = np.array([x0, y0, k], dtype=float)

    mu = 0.1
    gamma = 0.2
    nu = 0.1
    params = [mu, gamma, nu]

    Ns = 1000
    ds = 0.03
    tol = 1e-10
    max_newton = 30

    # tangente inicial (opción A): mu_s = 1
    DmuF = np.array([-1.0, -1.0])
    J = jacobian(U[0], U[1], U[2], params)
    dX = - np.linalg.solve(J, DmuF)      # u_s = -A^{-1} b
    dU = np.array([dX[0], dX[1], 1.0])
    dU = dU / np.linalg.norm(dU)

    # predictor inicial
    sign = 1  # or -1 to flip direction
    U_pred = U + sign * dU * ds
    ks = []
    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    for i in range(Ns):
        # --- corrector (Newton extendido sobre [F1,F2,H]) ---
        Ucorr = U_pred.copy()
        for it in range(max_newton):
            J_full = monster_jacobian(Ucorr[0], Ucorr[1], Ucorr[2], params, dU)   # usa el tangente previo (unitario)
            RHS = H_fun(Ucorr[0], Ucorr[1], Ucorr[2], params, U_pred, dU)         # residual [F1,F2,H]
            delta = np.linalg.solve(J_full, -RHS)                         # *** clave: -RHS ***
            Ucorr = Ucorr + delta
            if np.linalg.norm(delta) < tol * (1.0 + np.linalg.norm(Ucorr)):
                break

        # punto corregido
        U = Ucorr.copy()

        # monitoreo (norma de F en el punto corregido)
        F = function(U[0], U[1], U[2], params)
        Jstab = jacobian(U[0], U[1], U[2], params)
        lam = np.linalg.eigvals(Jstab)
        re = np.real(lam)
        size = 10
        # clasificar
        if np.all(re < 0):
            kind = "stable"  # sólido
        elif np.all(re > 0):
            kind = "unstable"  # vacío
        else:
            kind = "saddle"  # X

        # dibujar (k vs phi1 y k vs phi2) con marcador según 'kind'
        if kind == "saddle":
            ax.scatter(U[2], U[0], marker="x", c="b", s=size)
            ax.scatter(U[2], U[1], marker="x", c="r", s=size)
        elif kind == "stable":
            ax.scatter(U[2], U[0], marker="o", c="b", s=size)  # sólido
            ax.scatter(U[2], U[1], marker="o", c="r", s=size)
        else:  # unstable
            ax.scatter(U[2], U[0], marker="o", facecolors="w", edgecolors="b", s=size)  # vacío
            ax.scatter(U[2], U[1], marker="o", facecolors="w", edgecolors="r", s=size)
        dU_prev = dU.copy()   # guardar tangente viejo

        # --- recomputa tangente en el nuevo punto ---
        J = jacobian(U[0], U[1], U[2], params)
        DmuF = np.array([-1.0, -1.0])             # b = ∂F/∂mu
        dX = np.linalg.solve(J, DmuF) * (-1.0)   # u_s = -A^{-1} b, con mu_s=1
        dU = np.array([dX[0], dX[1], 1.0])
        dU = dU / np.linalg.norm(dU)

        # alinear orientación del tangente
        if np.dot(dU, dU_prev) < 0:
            dU = -dU

        # --- predictor siguiente ---
        U_pred = U + sign * dU * ds
        ks.append(U_pred[2])
    ks = np.array(ks)
    #plt.xlim(- 2 * np.pi, 2 * np.pi)
    plt.ylim(-np.pi, np.pi)
    yticks = [-np.pi, -np.pi/2, 0, np.pi/2, np.pi]
    ytick_labels = [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"]
    ax.tick_params(axis="y", direction="in", labelsize=15, left=True, right=True, labelleft=True, labelright=False)
    ax.tick_params(axis="x", direction="in", labelsize=15, top=True, bottom=True, labeltop=False, labelbottom=True)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ytick_labels)
    ax.set_xlabel("$\kappa$", fontsize=18)
    ax.set_ylabel("$\\theta_i$", fontsize=18)
    ax.grid(alpha=0.3)

    fig.subplots_adjust(left=0.2, right=0.9, bottom=0.15, top=0.9)
    plt.savefig("recuerdodemiprimeracontinuacionnumerica.png", dpi=300)