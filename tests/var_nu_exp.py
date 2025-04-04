from back_process import *

if __name__ == '__main__':
    a_ref = 11
    f_ref = 14.2
    gamma_ref = a_ref * f_ref ** 2
    f = np.arange(14.2, 15, 0.05)
    nus = 0.5 * ((f / 13.5) ** 2 - 1)
    print(nus)
    a = gamma_ref / f ** 2
    plt.scatter(f, a, c="k")
    plt.plot(f, a, color="k")
    plt.show()