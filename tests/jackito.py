from functions import *
from back_process import *

if __name__ == '__main__':
    a = np.loadtxt("C:/Users/Rafael Riveros/Documents/csv_jack.csv", delimiter=";", dtype=str)
    premios = ["Collar"]#["Pack molinillos FS","Cosmetiquero YvesSaintLaurent","Caja FS 1","Caja FS 2","Polera+totebag stitch","Libro ciencia oscuta (Gabriel Leon)","Espumante Xtra Brut (Portal del Alto)","Espumante Brut (Baron Lacroix)","Cabernet Sauvignion (Miguel Torres)","Collar de piedritas","Collar","Set de tacitas de cafe", "Caja de Alfajores"]
    N_premios = len(premios)
    for i in range(N_premios):
        i_rand = np.random.randint(0, 399)
        while a[i_rand, 1] == '':
            i_rand = np.random.randint(0, 399)
        j_rand = np.random.randint(0, len(premios))
        print(a[i_rand, 1] + " ha ganado " + premios[j_rand])
        premios.remove(premios[j_rand])
