import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def fit_plane(years, kms, prices, query_year, query_km):
    # Fit a plane to the data
    A = np.column_stack((years, kms, np.ones_like(years)))
    coeffs, _, _, _ = np.linalg.lstsq(A, prices, rcond=None)

    # Extract coefficients
    a, b, c = coeffs

    # Calculate price at the requested point
    estimated_price = a * query_year + b * query_km + c
    return estimated_price, coeffs

from functions import *
from back_process import *
from time_integrators import *

if __name__ == '__main__':
    # Data from the image (years, kilometers, price)
    data = [
        [146000, 2005, 2000000],
        [134000, 2005, 3200000],
        [132000, 2007, 3600000],
        [115000, 2007, 2700000],
        [98000, 2008, 3000000],
        [172000, 2005, 1950000],
        [151000, 2006, 1900000],
        [134000, 2005, 2500000],
        [120000, 2009, 3800000],
        [98000, 2008, 3000000],
        [151000, 2006, 1900000],
        [172000, 2005, 1950000],
        [157000, 2005, 2000000],
    ]

    # Separate into columns for processing
    kms, years, prices = zip(*data)

    # Convert to numpy arrays
    kms = np.array(kms)
    years = np.array(years)
    prices = np.array(prices)


    print(prices)
    print(kms)
    plt.scatter(kms, prices)
    plt.show()
    plt.close()

    plt.scatter(years, prices)
    plt.show()
    plt.close()

    query_year = 2006
    query_km = 85000

    # Compute plane fit and estimated price
    estimated_price, coeffs = fit_plane(years, kms, prices, query_year, query_km)
    print(f"Estimated price: {estimated_price}")

    fig = plt.figure(figsize=(5, 4))
    ax = fig.add_subplot(111, projection='3d')

    # Scatter plot of data points
    ax.scatter(years, kms, prices, color='red', label='Data Points')
    ax.scatter(2006, 85000, estimated_price, c="k")
    ax.set_xlabel('Year')
    ax.set_ylabel('Kilometers')
    ax.set_zlabel('Price')
    ax.set_title('Car Price Plane Fit')
    ax.legend()

    plt.show()