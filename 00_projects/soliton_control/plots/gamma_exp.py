from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    #### gamma = 0.17 ####
    nu_01o = np.loadtxt('D:/mnustes_science/experimental_data/soliton_control/gamma=0.17/bajando/toma_1/nus.txt', delimiter=',')
    center_01o = np.loadtxt('D:/mnustes_science/experimental_data/soliton_control/gamma=0.17/bajando/toma_1/centers.txt', delimiter=',')

    nu_01a = np.loadtxt('D:/mnustes_science/experimental_data/soliton_control/gamma=0.17/bajando/toma_2/nus.txt', delimiter=',')
    center_01a = np.loadtxt('D:/mnustes_science/experimental_data/soliton_control/gamma=0.17/bajando/toma_2/centers.txt', delimiter=',')

    nu_01b = np.loadtxt('D:/mnustes_science/experimental_data/soliton_control/gamma=0.17/bajando/toma_3/nus.txt', delimiter=',')
    center_01b = np.loadtxt('D:/mnustes_science/experimental_data/soliton_control/gamma=0.17/bajando/toma_3/centers.txt', delimiter=',')

    nu_01c = np.loadtxt('D:/mnustes_science/experimental_data/soliton_control/gamma=0.17/bajando/toma_4/nus.txt', delimiter=',')
    center_01c = np.loadtxt('D:/mnustes_science/experimental_data/soliton_control/gamma=0.17/bajando/toma_4/centers.txt', delimiter=',')

    #nu_01a = np.loadtxt('D:/mnustes_science/experimental_data/soliton_control/gamma=0.17/subiendo/toma_1/nus.txt', delimiter=',')
    #center_01a = np.loadtxt('D:/mnustes_science/experimental_data/soliton_control/gamma=0.17/subiendo/toma_1/centers.txt', delimiter=',')

    #### gamma = 0.18 ####
    nu_02a = np.loadtxt('D:/mnustes_science/experimental_data/soliton_control/gamma=0.18/bajando/toma_5/nus.txt', delimiter=',')
    center_02a = np.loadtxt('D:/mnustes_science/experimental_data/soliton_control/gamma=0.18/bajando/toma_5/centers.txt', delimiter=',')
    nu_02b = np.loadtxt('D:/mnustes_science/experimental_data/soliton_control/gamma=0.18/bajando/toma_6/nus.txt', delimiter=',')
    center_02b = np.loadtxt('D:/mnustes_science/experimental_data/soliton_control/gamma=0.18/bajando/toma_6/centers.txt', delimiter=',')
    nu_02c = np.loadtxt('D:/mnustes_science/experimental_data/soliton_control/gamma=0.18/bajando/toma_7/nus.txt', delimiter=',')
    center_02c = np.loadtxt('D:/mnustes_science/experimental_data/soliton_control/gamma=0.18/bajando/toma_7/centers.txt', delimiter=',')

    #### gamma = 0.19 ####
    nu_03a = np.loadtxt('D:/mnustes_science/experimental_data/soliton_control/gamma=0.19/subiendo/nus.txt', delimiter=',')
    center_03a = np.loadtxt('D:/mnustes_science/experimental_data/soliton_control/gamma=0.19/subiendo/centers.txt', delimiter=',')
    nu_03b = np.loadtxt('D:/mnustes_science/experimental_data/soliton_control/gamma=0.19/bajando/toma_1/nus.txt', delimiter=',')
    center_03b = np.loadtxt('D:/mnustes_science/experimental_data/soliton_control/gamma=0.19/bajando/toma_1/centers.txt', delimiter=',')
    nu_03c = np.loadtxt('D:/mnustes_science/experimental_data/soliton_control/gamma=0.19/bajando/toma_2/nus.txt', delimiter=',')
    center_03c = np.loadtxt('D:/mnustes_science/experimental_data/soliton_control/gamma=0.19/bajando/toma_2/centers.txt', delimiter=',')
    nu_03d = np.loadtxt('D:/mnustes_science/experimental_data/soliton_control/gamma=0.19/bajando/toma_3/nus.txt', delimiter=',')
    center_03d = np.loadtxt('D:/mnustes_science/experimental_data/soliton_control/gamma=0.19/bajando/toma_3/centers.txt', delimiter=',')
    nu_03e = np.loadtxt('D:/mnustes_science/experimental_data/soliton_control/gamma=0.19/bajando/toma_4/nus.txt', delimiter=',')
    center_03e = np.loadtxt('D:/mnustes_science/experimental_data/soliton_control/gamma=0.19/bajando/toma_4/centers.txt', delimiter=',')

    center_02a_mask = []
    for i in range(len(nu_02b)):
        if i % 2 == 0:
            center_02a_mask.append(center_02a[int(i / 2)])
        else:
            center_02a_mask.append(np.nan)
    center_02 = np.nanmean(np.array([center_02a_mask, center_02b, center_02c]), axis=0)
    center_std_02 = np.nanstd(np.array([center_02a_mask, center_02b, center_02c]), axis=0)

    center_03 = np.nanmean(np.array([center_03b, center_03c]), axis=0)
    center_std_03 = np.nanstd(np.array([center_03b, center_03c]), axis=0)

    #center_01 = np.nanmean(np.array([center_01a, center_01b, center_01c, center_01o]), axis=0)
    #center_std_01 = np.nanstd(np.array([center_01a, center_01b, center_01c, center_01o]), axis=0)

    #center_01a = np.append(center_01a, np.zeros(len(np.arange(nu_02b[0] + 0.005, 0, 0.005))))
    #nu_01a = np.append(nu_01a, np.arange(nu_01a[0] + 0.005, 0, 0.005))

    #center_02 = np.append(np.zeros(len(np.arange(nu_02b[0] + 0.005, 0, 0.005))), center_02)
    #center_std_02 = np.append(np.zeros(len(np.arange(nu_02b[0] + 0.005, 0, 0.005))), center_std_02)
    #nu_02b = np.append(np.flip(np.arange(nu_02b[0] + 0.005, 0, 0.005)), nu_02b)
    #plt.scatter(nu_01o, np.abs(center_01o), c="y", label="01")# label="$\gamma=0.17$")
    #plt.scatter(nu_01a, np.abs(center_01a), c="g", label="02")# label="$\gamma=0.17$")
    #plt.scatter(nu_01b, np.abs(center_01b), c="r", label="03")# label="$\gamma=0.17$")
    #plt.scatter(nu_01c, np.abs(center_01c), c="b", label="04")# label="$\gamma=0.17$")
    plt.scatter(nu_03a, np.abs(center_03a), c="g", label="01")# label="$\gamma=0.17$")
    plt.scatter(nu_03b, np.abs(center_03b), c="k", label="02")# label="$\gamma=0.17$")
    plt.scatter(nu_03d, np.abs(center_03d), c="r", label="03")# label="$\gamma=0.17$")
    plt.scatter(nu_03e, np.abs(center_03e), c="b", label="04")# label="$\gamma=0.17$")
    #plt.scatter(nu_03a, np.abs(center_03a))
    #plt.scatter(nu_03a, np.abs(center_03a))
    #plt.errorbar(nu_02b, np.abs(center_02), center_std_02, marker="o", ms=6, label="$\gamma=0.18$")
    #plt.errorbar(nu_03a, np.abs(center_03), center_std_03, marker="o", ms=6, label="$\gamma=0.19$")
    #plt.scatter(nu_03c, np.abs(center_03c), c="gold", label="$\gamma=0.19$")
    plt.legend()
    plt.show()