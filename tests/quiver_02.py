import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    # Creating arrow
    x_min = -1
    x_max = - x_min
    y_min = x_min
    y_max = x_max
    x = [0] #np.arange(x_min, x_max, 0.05)
    y = np.arange(y_min, y_max, 0.05)
    z = np.arange(y_min, y_max, 0.05)

    X, Y, Z = np.meshgrid(x, y, z)

    X_dot = - X - Y - Z - Y * Z
    Y_dot = - Z - Y - Z - X * Z
    Z_dot = - X - Y - Z - X * Y
    # creating plot
    ax = plt.figure().add_subplot(projection='3d')
    ax.quiver(X, Y, Z, X_dot, Y_dot, X_dot, length=0.05, normalize=True, arrow_length_ratio=0.3)

    plt.grid(alpha=0.3)
    ax.set_aspect('equal')

    # show plot
    plt.show()