import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D, art3d


if __name__ == '__main__':
    # Create random data
    np.random.seed(0)
    x = np.random.normal(0, 1, 100)
    y = np.random.normal(0, 1, 100)
    z = np.random.normal(0, 1, 100)
    colors = np.sqrt(x ** 2 + y ** 2 + z ** 2)

    # Set up plot
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    # Scatter plot
    scatter = ax.scatter(
        x, y, z,
        c=colors,
        cmap='plasma',
        s=60,
        alpha=1.0,
        edgecolors='white'
    )

    # Label axes
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # Remove pane colors (ugly box)
    for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
        axis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        #axis._axinfo['grid'].update({'linewidth': 0})  # Remove grid

    # Get axis limits
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    zlim = ax.get_zlim()


    # Draw black lines for X, Y, Z axes manually (outer frame)
    def draw_line(p1, p2):
        line = art3d.Line3D(*zip(p1, p2), color='black', linewidth=1.5)
        ax.add_line(line)


    # Bottom square
    draw_line((xlim[0], ylim[0], zlim[0]), (xlim[1], ylim[0], zlim[0]))
    draw_line((xlim[0], ylim[0], zlim[0]), (xlim[0], ylim[1], zlim[0]))
    draw_line((xlim[1], ylim[0], zlim[0]), (xlim[1], ylim[1], zlim[0]))
    draw_line((xlim[0], ylim[1], zlim[0]), (xlim[1], ylim[1], zlim[0]))

    # Vertical edges
    draw_line((xlim[0], ylim[0], zlim[0]), (xlim[0], ylim[0], zlim[1]))
    draw_line((xlim[1], ylim[0], zlim[0]), (xlim[1], ylim[0], zlim[1]))
    draw_line((xlim[1], ylim[1], zlim[0]), (xlim[1], ylim[1], zlim[1]))
    #draw_line((xlim[0], ylim[1], zlim[0]), (xlim[0], ylim[1], zlim[1]))

    # Top square
    draw_line((xlim[0], ylim[0], zlim[1]), (xlim[1], ylim[0], zlim[1]))
    #draw_line((xlim[0], ylim[0], zlim[1]), (xlim[0], ylim[1], zlim[1]))
    draw_line((xlim[1], ylim[1], zlim[1]), (xlim[1], ylim[0], zlim[1]))
    #draw_line((xlim[1], ylim[1], zlim[1]), (xlim[0], ylim[1], zlim[1]))

    # Set view
    ax.view_init(elev=30, azim=135)
    ax.grid()
    # Colorbar
    cbar = fig.colorbar(scatter, shrink=0.6)
    cbar.set_label('Distance from origin')

    plt.title('Pretty 3D Scatter Plot with Black Axes', fontsize=14)
    plt.tight_layout()
    plt.show()