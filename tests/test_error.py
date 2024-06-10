from time_integrators import *

if __name__ == '__main__':
    F = [1]
    n = [1]
    phase = 0.1
    for k in n:
        s_k = np.sin(2 * np.pi * k * phase)
        c_k = np.cos(2 * np.pi * k * phase)
        F.append(s_k)
        F.append(c_k)
    A = np.outer(F, F)
    print(F)
    print(A)
    det_A = np.linalg.det(A)
    print(det_A)