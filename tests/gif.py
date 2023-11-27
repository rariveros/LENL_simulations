from back_process import *
from matplotlib.animation import FuncAnimation, PillowWriter

if __name__ == "__main__":

    ### Definiendo parametros y eligiendo carpeta a detectar ###
    disco = 'D:'
    project_file = 'OSC'
    initial_dir_img = str(disco) + '/mnustes_science/images'
    initial_dir_data = str(disco) + '/mnustes_science/experimental_data'
    root = tk.Tk()
    root.withdraw()
    file = filedialog.askdirectory(parent=root, initialdir=initial_dir_img, title='Elección de carpeta')
    parent_file_name = os.path.basename(file)
    IMG_names = os.listdir(file)
    N_img = len(IMG_names)
    imgs_strobo = []
    for i_strobo in range(N_img):
        img_i = cv2.imread(file + '/' + IMG_names[int(i_strobo)])
        #image_rgb_i = cv2.cvtColor(img_i, cv2.COLOR_BGR2RGB)
        imgs_strobo.append(img_i)

    fig = plt.figure()
    im = plt.imshow(imgs_strobo[0], animated=True)
    plt.axis('off')

    def updatefig(i):
        im.set_array(imgs_strobo[i])
        #plt.title("$t = " + str(T_strobo[i]).split(".")[0] + "." + str(T_strobo[i]).split(".")[1][0:1] + "\ \\textrm{s}$", fontsize=25)
        return im,



    ani = FuncAnimation(fig, updatefig, frames=N_img, interval=1)
    FFwriter = animation.FFMpegWriter()
    ani.save('name_test.gif', writer=FFwriter, dpi=300)
    plt.close()