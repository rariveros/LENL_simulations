from functions import *
from back_process import *
from time_integrators import *
from scipy.sparse import diags
from scipy.sparse.linalg import splu

if __name__ == '__main__':
    def sparse_DD_noborder(Nx, dx):
        data = np.ones((3, Nx))
        data[1] = -2 * data[1]
        diags = [-1, 0, 1]
        D2 = sparse.spdiags(data, diags, Nx, Nx) / (dx ** 2)
        D2 = sparse.lil_matrix(D2)
        return D2
    dt = 1
    T = 100
    [tmin, tmax, dt] = [0, T, dt]
    t_grid = np.arange(tmin, tmax + dt, dt)
    freq = 1
    epsilon = 0.00001
    Nt = len(t_grid)
    D2 = sparse_DD_noborder(Nt, dt)
    data = np.ones((1, Nt))
    data[0] =  freq ** 2 * data[0]
    diags = [0]
    FREQ = sparse.spdiags(data, diags, Nt, Nt)
    FREQ = sparse.lil_matrix(FREQ)
    matrix_HO = -(D2 + epsilon ** 2 * FREQ)
    matrix_FREE = -(D2)
    det = np.linalg.det(matrix_HO.todense())
    det_0 = np.linalg.det(matrix_FREE.todense())
    print(det)
    freq_aux = (2 / epsilon) * np.arcsin(epsilon * freq / 2)
    xd = np.sin(freq * T) / freq
    print(xd)