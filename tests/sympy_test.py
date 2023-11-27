from sympy import *

if __name__ == '__main__':
    # Creating arrow
    u = Symbol('u')
    v = Symbol('v')
    w = Symbol('w')

    a1 = Symbol('a1')
    b1 = Symbol('b1')
    c1 = Symbol('c1')
    a2 = Symbol('a2')
    b2 = Symbol('b2')
    c2 = Symbol('c2')

    h1 =a1 * u ** 2 + b1 * u ** 3 + a2 * v ** 2 + b2 * v ** 3
    Dh1 = diff(h1, w)
    h2 = a2 * w ** 2 + b2 * w ** 3# + c2 * w ** 4#c2 * w ** 4 #
    Dh2 = diff(h2, w)
    def F(u, v, w):
        x = (- u - v + w)
        y = (v + w)
        z = (u + w)
        f1 = (-y * z - x * z + 2 * y * z) / 3
        f2 = (-y * z + 2 * x * z - y * z) / 3
        f3 = (y * z + x * z + y * z) / 3
        return f1, f2, f3

    A = 0
    B = C = 3
    f1, f2, f3 = F(u, v, h1)
    N = (Dh1) * (A * u + A * v + f1 + f2) - B * h1 - f3
    V = [expand(N).coeff(w**1), expand(N).coeff(w**2), expand(N).coeff(w**3), expand(N).coeff(w**4), expand(N).coeff(w**5)]
    print(solve(V))

    g1, g2, g3 = F(((4/27) * w ** 2 - (2 / 243) * w ** 3), -(4/27) * w ** 2 - (46 / 243) * w ** 3, w)
    X = - 3 * u + g1
    Y = - 3 * v + g2
    Z = g3
    print(expand(Z))

