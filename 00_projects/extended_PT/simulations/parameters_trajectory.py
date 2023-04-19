from time_integrators import *

if __name__ == '__main__':

    pist_dist = np.arange(4, 7.5, 0.5)
    f_i = np.arange(15, 17.25, 0.25)
    a = np.arange(10, 13)
    params = []
    params_exp = []
    for i in range(len(f_i)):
        params_i = []
        params_exp_i = []
        for j in range(len(a)):
            alpha, beta, nu, gamma_0 = fluid_pdnls_parameters(f_i[i], a[j], d=20)
            print(f_i[i])
            print(nu)
            params_i.append([alpha, beta, nu, gamma_0])
            params_exp_i.append([f_i[i], a[j]])
        params.append(params_i)
        params_exp.append(params_exp_i)
    params = np.array(params)
    params_exp = np.array(params_exp)


    for k in range(len(f_i)):
        plt.plot(params[k, :, 2], params[k, :, 3], color="k")
        plt.scatter(params[k, :, 2], params[k, :, 3], c="k")
    plt.show()
    plt.close()

    for k in range(len(f_i)):
        plt.plot(params_exp[k, :, 0], params_exp[k, :, 1], color="k")
        plt.scatter(params_exp[k, :, 0], params_exp[k, :, 1], c="k")
    plt.show()