import matplotlib.pyplot as plt

from functions import *
from back_process import *
from time_integrators import *

import numpy as np


def generate_discrete_function(x_grid, left_square_width, left_square_amplitude, right_square_width,
                               right_square_amplitude):
    f_i = np.zeros_like(x_grid)

    # Calculate separation to center it in the middle of x_grid
    separation = (x_grid[-1] - x_grid[0] - left_square_width - right_square_width) / 2

    # Define left square
    left_square_start = x_grid[int(len(x_grid) / 4)]
    left_square_end = left_square_start + left_square_width
    f_i[(x_grid >= left_square_start) & (x_grid <= left_square_end)] = left_square_amplitude

    # Define right square
    right_square_start = x_grid[-1] - right_square_width - separation
    right_square_end = right_square_start + right_square_width
    f_i[(x_grid >= right_square_start) & (x_grid <= right_square_end)] = right_square_amplitude

    return f_i


if __name__ == '__main__':
    # Example usage:
    n_points = 100
    x_grid = np.linspace(0, 100, n_points)
    left_square_width = 10
    left_square_amplitude = 5
    right_square_width = 15
    right_square_amplitude = 8

    f_i = generate_discrete_function(x_grid, left_square_width, left_square_amplitude, right_square_width,
                                     right_square_amplitude)
    plt.plot(x_grid, f_i)
    plt.show()