from back_process import *
from itertools import zip_longest


if __name__ == "__main__":
    test_directory = "D:/mnustes_science/simulation_data/FD/aa"

    # HACER SIMULACION LARGA DE ESTE Y TOMAR TODOS LOS PATRONES DESDE EL MISMO T, LUEGO BINEAR EN X
    mode = "solitons" #solitons or fronts

    root = tk.Tk()
    root.withdraw()
    directory = filedialog.askdirectory(parent=root, initialdir=test_directory, title='Elección de carpeta')
    parent_directory_name = os.path.basename(directory)
    save_directory = directory + "//analysis"

    root = tk.Tk()
    root.withdraw()
    reference_image = filedialog.askopenfilename(parent=root, initialdir=directory, title='Reference Selection')
    img_reference = cv2.imread(str(reference_image))
    data_type = "SIM"
    if data_type == "SIM":
        Z_r = np.loadtxt(directory + '/field_real.txt', delimiter=',')
        Z_i = np.loadtxt(directory + '/field_img.txt', delimiter=',')
        x_grid = np.loadtxt(directory + '/X.txt', delimiter=',')
        y_grid = np.loadtxt(directory + '/T.txt', delimiter=',')
        Z_complex = Z_r + 1j * Z_i
        Z = np.absolute(Z_complex)
        x_grid = x_grid.astype('float64')
        Z = Z.astype('float64')
    elif data_type == "EXP":
        Z = np.loadtxt(directory + '/Z_strobo.txt', delimiter=',')
        x_grid = np.loadtxt(directory + '/X_mm.txt', delimiter=',')
        y_grid = np.loadtxt(directory + '/T_strobo.txt', delimiter=',')
    elif data_type == "ELRAM":
        A = np.loadtxt('F:/SpatioTemporalX.dat', delimiter='\t')
        X = A[:, 0]
        T = A[:, 1]
        Z = A[:, 2]
        Nx = 110
        Nt = 9999
        Z = Z[1200:].reshape((Nt, Nx))
        x_grid = X[1200:1200 + 110]
        y_grid = T[1200::110]
    Nx = len(x_grid)
    Ny = len(y_grid)
    dx = x_grid[1] - x_grid[0]
    dy = y_grid[1] - y_grid[0]

    resize_scale = 0.4
    h, w, c = img_reference.shape
    h_resized, w_resized = h * resize_scale, w * resize_scale
    resized_img = cv2.resize(img_reference, (int(w_resized), int(h_resized)))
    cut_coords = cv2.selectROI(resized_img)
    cv2.destroyAllWindows()

    ### Cut image with resized scale ###
    cut_coords_list = list(cut_coords)
    x_1 = int(cut_coords_list[0] / resize_scale)
    x_2 = int(cut_coords_list[2] / resize_scale)
    y_1 = int(cut_coords_list[1] / resize_scale)
    y_2 = int(cut_coords_list[3] / resize_scale)
    img_crop = img_reference[y_1:(y_1 + y_2), x_1:(x_1 + x_2)]
    img_gray_ref = cv2.cvtColor(img_crop, cv2.COLOR_BGR2GRAY)

    Ny_pix, Nx_pix = img_gray_ref.shape
    N_structures = 1
    window_info = []
    for k in range(N_structures):
        wind_R_pix, wind_L_pix, t_init_pix = define_window(img_crop, resize_scale)
        wind_L_j = int((wind_L_pix / Nx_pix) * Nx)
        wind_R_j = int((wind_R_pix / Nx_pix) * Nx)
        t_init_i = int(Ny - (t_init_pix / Ny_pix) * Ny)
        delta_wind_0 = wind_R_j - wind_L_j
        if delta_wind_0 % 2 != 0:
            wind_R_j = wind_R_j + 1
        delta_wind_0 = wind_R_j - wind_L_j
        window_info.append([delta_wind_0, wind_L_j, wind_R_j, t_init_i])
    Js = []
    Is = []
    Xs = []
    Ys = []
    Zs = []
    # HAY QUE AGREGAR UN EXIT EN EL FOR DEL I SI LA ESTRUCTURA DESAPARECE
    for k in range(N_structures):
        [delta_wind_0, wind_L_j, wind_R_j, t_init_i] = window_info[k]
        Z_windowed_0 = Z[0, wind_L_j:wind_R_j]
        N_window = len(Z_windowed_0)
        D = sparse_D_neumann_4order(N_window, dx)
        DD = sparse_DD(N_window, dx)
        J = []
        I = []
        X_structure = []
        Y_structure = []
        Z_structure = []
        i_init = int(t_init_i)
        X_max_i = 0
        for i in range(i_init, Ny):
            if i != i_init:
                wind_L_j = wind_L_j_new
                wind_R_j = wind_R_j_new
            Z_windowed_i = Z[i, wind_L_j:wind_R_j]
            if mode == "solitons":
                j_max_i = np.argmax(Z_windowed_i)
            elif mode == "fronts":
                DZ = Dx(D, Z_windowed_i)
                j_max_i = np.argmax(DZ)
            if j_max_i == 0 or j_max_i == len(Z_windowed_i) or j_max_i == len(Z_windowed_i) + 1 or np.isnan(X_max_i):
                X_max_i = float("nan")
                Z_max_i = float("nan")
            else:
                ajuste = np.polyfit(dx * x_grid[wind_L_j:wind_R_j], Z[i, wind_L_j:wind_R_j], 2)
                X_max_i = - 2 * ajuste[1] / ajuste[0]
                #X_max_i = x_grid[wind_L_j + j_max_i]
                Z_max_i = Z[i, wind_L_j + j_max_i]
            Y_max_i = y_grid[i]
            J.append(j_max_i)
            I.append(i)
            X_structure.append(X_max_i)
            Y_structure.append(Y_max_i)
            Z_structure.append(Z_max_i)
            wind_L_j_new = int(wind_L_j + j_max_i - int(delta_wind_0 / 2))
            wind_R_j_new = int(wind_L_j + j_max_i + int(delta_wind_0 / 2))
        J = np.array(J)
        I = np.array(I)
        X_structure = np.array(X_structure)
        Y_structure = np.array(Y_structure)
        Z_structure = np.array(Z_structure)
        Js.append(J)
        Is.append(I)
        Xs.append(X_structure)
        Ys.append(Y_structure)
        Zs.append(Z_structure)
    Js = np.array([list(tpl) for tpl in zip(*zip_longest(*Js))], dtype=np.float64)
    Is = np.array([list(tpl) for tpl in zip(*zip_longest(*Is))], dtype=np.float64)
    Xs = np.array([list(tpl) for tpl in zip(*zip_longest(*Xs))], dtype=np.float64)
    Ys = np.array([list(tpl) for tpl in zip(*zip_longest(*Ys))], dtype=np.float64)
    Zs = np.array([list(tpl) for tpl in zip(*zip_longest(*Zs))], dtype=np.float64)
    Ny_effective = len(Ys[0])
    Xs = filtro_superficie(Xs, 200, "X")
    Dy = sparse_D_neumann_4order(Ny_effective, dy)
    DX = Dx(Dy, np.transpose(Xs))
    mask_1 = (DX > 0.125)
    mask_2 = (DX < 0.00)
    DX[mask_1] = np.nan
    DX[mask_2] = np.nan

    delta_x_L = np.abs(Xs - (-37))
    delta_x_R = np.abs(Xs -(37))
    gaussian_L_j = np.nanargmin(delta_x_L[:])
    gaussian_R_j = np.nanargmin(delta_x_R[:])

    print(Xs[0, gaussian_L_j])
    print(Xs[0, gaussian_R_j])
    print(DX[gaussian_L_j][0])
    print(DX[gaussian_R_j][0])

    if not os.path.exists(save_directory):
        os.makedirs(save_directory)
    np.savetxt(save_directory + '//J.txt', Js, delimiter=',')
    np.savetxt(save_directory + '//I.txt', Is, delimiter=',')
    np.savetxt(save_directory + '//X_structure.txt', Xs, delimiter=',')
    np.savetxt(save_directory + '//Y_structure.txt', Ys, delimiter=',')
    np.savetxt(save_directory + '//Z_structure.txt', Zs, delimiter=',')

    ### Visualizacion del diagrama espacio-temporal  ###
    pcm = plt.pcolormesh(x_grid, y_grid, Z, cmap=parula_map, shading='auto')
    cbar = plt.colorbar(pcm, shrink=1)
    cbar.set_label('$|A(x, t)|$', rotation=0, size=20, labelpad=-27, y=1.1)
    plt.xlabel('$x$', size='20')
    plt.ylabel('$t$', size='20')
    plt.grid(linestyle='--', alpha=0.5)
    plt.plot(np.transpose(Xs), np.transpose(Ys), color="r")
    plt.savefig(save_directory + '//field_tracked.png', dpi=200)
    plt.close()

    plt.plot(np.transpose(Ys), np.transpose(Xs), color="k")
    plt.xlabel('$t\ (\\textrm{s})$', size='20')
    plt.ylabel('$x\ (\\textrm{mm})$', size='20')
    plt.grid(linestyle='--', alpha=0.5)
    plt.savefig(save_directory + '//crest_position.png', dpi=300)
    plt.close()

    plt.plot(np.transpose(Ys), DX, color="k")
    plt.xlabel('$t\ (\\textrm{s)}$', size='20')
    plt.ylabel('$v\ (\\textrm{mm/s})$', size='20')
    #plt.ylim(0, 0.2)
    plt.grid(linestyle='--', alpha=0.5)
    plt.savefig(save_directory + '//crest_velocity.png', dpi=300)
    plt.close()

    plt.plot(np.transpose(Ys), np.transpose(Zs), color="k")
    plt.xlabel('$t\ (\\textrm{s})$', size='20')
    plt.ylabel('$|A|\ (\\textrm{mm})$', size='20')
    plt.grid(linestyle='--', alpha=0.5)
    plt.savefig(save_directory + '//ampd_tracked.png', dpi=300)
    plt.close()

    plt.plot(np.transpose(Xs), np.transpose(Zs), color="k")
    plt.xlabel('$x\ (\\textrm{mm})$', size='20')
    plt.ylabel('$|A|\ (\\textrm{mm})$', size='20')
    plt.grid(linestyle='--', alpha=0.5)
    plt.savefig(save_directory + '//ampd(x).png', dpi=300)
    plt.close()

    plt.scatter(np.transpose(Xs), DX, color="k", s=1)
    plt.xlabel('$x\ (\\textrm{mm)}$', size='20')
    plt.ylabel('$v\ (\\textrm{mm/s})$', size='20')
    #plt.ylim(0.05, 0.125)
    plt.grid(linestyle='--', alpha=0.5)
    plt.savefig(save_directory + '//v(x).png', dpi=300)
    plt.close()