from back_process import *
from jacobians import *

def gaussian(x, A, x_0, sigma):
    return A * np.exp(- (x - x_0) ** 2 / (2 * sigma ** 2))

if __name__ == '__main__':
    eq = 'pdnlS'
    disco = 'D:/'
    #dists = ["20.000", "22.000", "24.000", "26.000", "28.000"]
    dists = ["30.2930", "30.2950", "30.2970", "30.3000", "30.3200", "30.3400", "30.3600", "30.3800",  "30.4000", "30.5000",  "30.7000", "30.8000", "30.9000", "31.000", "32.000",
             "33.000", "34.000","35.000", "36.000","37.000", "38.000","39.000", "40.000","41.000", "42.000","43.000", "44.000","45.000", "46.000","47.000", "48.000","49.000", "50.000","51.000", "52.000", "53.000","54.000","55.000", "56.000","57.000", "58.000","59.000"]
    #dists = ["32.000", "34.000", "36.000", "38.000", "40.000", "42.000", "44.000", "46.000", "48.000", "50.000", "52.000",
    #         "54.000", "56.000", "58.000"]
    colors = np.arange(0, 1, 1 / (len(dists) + 1))
    print(colors)
    real_maxs = []
    imag_maxs = []
    DIST = []
    parent_directory ="C:/mnustes_science/simulation_data/FD/PT_extended/variando_distancia/alpha=6.5240/beta=1.000/mu=0.1000/nu=0.0180/sigma=6.000/comprimido"
    dist_dir = [name for name in os.listdir(parent_directory) if os.path.isdir(os.path.join(parent_directory, name))]
    for i in range(len(dist_dir)):
        print(" #################   " + dist_dir[i] + "    ###########")
        directory = parent_directory + "/" + dist_dir[i]
        os.remove(directory + '/field_real.txt')
        os.remove(directory + '/field_img.txt')
