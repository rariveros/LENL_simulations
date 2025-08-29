from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    def function(x, y, mu):
        F1 = x ** 2 - mu
        F2 = y ** 2 - 1
        return np.array([F1, F2])

    def jacobian(x, y, mu):
        J11 = 2 * x
        J12 = 0
        J21 = 0
        J22 = 2 * y
        return np.array([[J11, J12], [J21, J22]])

    def H_fun(x, y, mu, U_pred, dU_s):
        x_pred = U_pred[0]
        y_pred = U_pred[1]
        mu_pred = U_pred[2]

        dx_s = dU_s[0]
        dy_s = dU_s[1]
        dm_s = dU_s[2]

        F1 = x ** 2 - mu
        F2 = y ** 2 - 1
        H = (x - x_pred) * dx_s + (y - y_pred) * dy_s + (mu - mu_pred) * dm_s
        return np.array([F1, F2, H])

    def monster_jacobian(x, y, mu, dU_s):
        dx_s = dU_s[0]
        dy_s = dU_s[1]
        dm_s = dU_s[2]
        J11 = 2 * x
        J12 = 0
        J13 = -1
        J21 = 0
        J22 = 2 * y
        J23 = 0
        J31 = dx_s
        J32 = dy_s
        J33 = dm_s
        return np.array([[J11, J12, J13], [J21, J22, J23], [J31, J32, J33]])

    # --- semilla ---
    mu = 1.0
    x0 = -np.sqrt(mu)+0.08
    y0 = 1.06
    U = np.array([x0, y0, mu], dtype=float)

    Ns = 1000
    ds = 0.01
    tol = 1e-10
    max_newton = 20

    # tangente inicial (opción A): mu_s = 1
    DmuF = np.array([-1.0, 0.0])
    J = jacobian(U[0], U[1], U[2])
    dX = - np.linalg.solve(J, DmuF)      # u_s = -A^{-1} b
    dU = np.array([dX[0], dX[1], 1.0])
    dU = dU / np.linalg.norm(dU)

    # predictor inicial
    sign = -1  # or -1 to flip direction
    U_pred = U + sign * dU * ds
    mus = []

    for i in range(Ns):
        # --- corrector (Newton extendido sobre [F1,F2,H]) ---
        Ucorr = U_pred.copy()
        for it in range(max_newton):
            J_full = monster_jacobian(Ucorr[0], Ucorr[1], Ucorr[2], dU)   # usa el tangente previo (unitario)
            RHS = H_fun(Ucorr[0], Ucorr[1], Ucorr[2], U_pred, dU)         # residual [F1,F2,H]
            delta = np.linalg.solve(J_full, -RHS)                         # *** clave: -RHS ***
            Ucorr = Ucorr + delta
            if np.linalg.norm(delta) < tol * (1.0 + np.linalg.norm(Ucorr)):
                break

        # punto corregido
        U = Ucorr.copy()

        # monitoreo (norma de F en el punto corregido)
        F = function(U[0], U[1], U[2])
        plt.scatter(U[2], U[0], c="b")
        plt.scatter(U[2], U[1], c="r")
        dU_prev = dU.copy()   # guardar tangente viejo
        # --- recomputa tangente en el nuevo punto ---
        J = jacobian(U[0], U[1], U[2])

        # --- recomputa tangente en el nuevo punto ---
        J = jacobian(U[0], U[1], U[2])
        DmuF = np.array([-1.0, 0.0])             # b = ∂F/∂mu
        dX = np.linalg.solve(J, DmuF) * (-1.0)   # u_s = -A^{-1} b, con mu_s=1
        dU = np.array([dX[0], dX[1], 1.0])
        dU = dU / np.linalg.norm(dU)


        # alinear orientación del tangente
        if np.dot(dU, dU_prev) < 0:
            dU = -dU

        # --- predictor siguiente ---
        U_pred = U + sign * dU * ds
        mus.append(U_pred[2])
    mus = np.array(mus)
    plt.plot(mus, np.sqrt(mus))
    plt.show()