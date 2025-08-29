import numpy as np

if __name__ == '__main__':
    # blanco → verde
    white = np.array([1.0, 1.0, 1.0])
    green = np.array([91/255, 249/255, 149/255])  # #5bf995

    N = 100
    gamma = 0.8  # más detalle cerca de 0
    x = np.linspace(0, 1, N) ** gamma

    custom_green_rgb_100 = [
        tuple(((1-xi) * white + xi * green).round(6)) for xi in x
    ]

    # ejemplo
    print(custom_green_rgb_100)
