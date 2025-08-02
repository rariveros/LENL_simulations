#### MATPLOTLIB ###
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
#import pandas as pd
from matplotlib import rc
from matplotlib.animation import FuncAnimation
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import FuncFormatter, MultipleLocator, FormatStrFormatter
from matplotlib.colors import TwoSlopeNorm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

import warnings
#warnings.filterwarnings("ignore", category=SyntaxWarning)

#plt.rcParams['animation.ffmpeg_path'] = "C:/Users/rariv/ffmpeg/bin/ffmpeg.exe"
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

import cv2

#### SCIPY ###
from scipy.signal import filtfilt
from scipy.optimize import curve_fit
from scipy.integrate import odeint
import scipy.sparse as sparse
from scipy import special
from scipy import signal
from scipy.fft import fft, fftfreq
from scipy.signal import hilbert, chirp
from scipy.integrate import solve_ivp
from scipy.sparse import diags
from scipy import integrate
from scipy.optimize import fsolve

#### OTROS ###
import numpy as np
import os
import time
import datetime
import shutil
import tkinter as tk
from tkinter import filedialog


cm_data = [[0.2422, 0.1504, 0.6603],
[0.2444, 0.1534, 0.6728],
[0.2464, 0.1569, 0.6847],
[0.2484, 0.1607, 0.6961],
[0.2503, 0.1648, 0.7071],
[0.2522, 0.1689, 0.7179],
[0.254, 0.1732, 0.7286],
[0.2558, 0.1773, 0.7393],
[0.2576, 0.1814, 0.7501],
[0.2594, 0.1854, 0.761],
[0.2611, 0.1893, 0.7719],
[0.2628, 0.1932, 0.7828],
[0.2645, 0.1972, 0.7937],
[0.2661, 0.2011, 0.8043],
[0.2676, 0.2052, 0.8148],
[0.2691, 0.2094, 0.8249],
[0.2704, 0.2138, 0.8346],
[0.2717, 0.2184, 0.8439],
[0.2729, 0.2231, 0.8528],
[0.274, 0.228, 0.8612],
[0.2749, 0.233, 0.8692],
[0.2758, 0.2382, 0.8767],
[0.2766, 0.2435, 0.884],
[0.2774, 0.2489, 0.8908],
[0.2781, 0.2543, 0.8973],
[0.2788, 0.2598, 0.9035],
[0.2794, 0.2653, 0.9094],
[0.2798, 0.2708, 0.915],
[0.2802, 0.2764, 0.9204],
[0.2806, 0.2819, 0.9255],
[0.2809, 0.2875, 0.9305],
[0.2811, 0.293, 0.9352],
[0.2813, 0.2985, 0.9397],
[0.2814, 0.304, 0.9441],
[0.2814, 0.3095, 0.9483],
[0.2813, 0.315, 0.9524],
[0.2811, 0.3204, 0.9563],
[0.2809, 0.3259, 0.96],
[0.2807, 0.3313, 0.9636],
[0.2803, 0.3367, 0.967],
[0.2798, 0.3421, 0.9702],
[0.2791, 0.3475, 0.9733],
[0.2784, 0.3529, 0.9763],
[0.2776, 0.3583, 0.9791],
[0.2766, 0.3638, 0.9817],
[0.2754, 0.3693, 0.984],
[0.2741, 0.3748, 0.9862],
[0.2726, 0.3804, 0.9881],
[0.271, 0.386, 0.9898],
[0.2691, 0.3916, 0.9912],
[0.267, 0.3973, 0.9924],
[0.2647, 0.403, 0.9935],
[0.2621, 0.4088, 0.9946],
[0.2591, 0.4145, 0.9955],
[0.2556, 0.4203, 0.9965],
[0.2517, 0.4261, 0.9974],
[0.2473, 0.4319, 0.9983],
[0.2424, 0.4378, 0.9991],
[0.2369, 0.4437, 0.9996],
[0.2311, 0.4497, 0.9995],
[0.225, 0.4559, 0.9985],
[0.2189, 0.462, 0.9968],
[0.2128, 0.4682, 0.9948],
[0.2066, 0.4743, 0.9926],
[0.2006, 0.4803, 0.9906],
[0.195, 0.4861, 0.9887],
[0.1903, 0.4919, 0.9867],
[0.1869, 0.4975, 0.9844],
[0.1847, 0.503, 0.9819],
[0.1831, 0.5084, 0.9793],
[0.1818, 0.5138, 0.9766],
[0.1806, 0.5191, 0.9738],
[0.1795, 0.5244, 0.9709],
[0.1785, 0.5296, 0.9677],
[0.1778, 0.5349, 0.9641],
[0.1773, 0.5401, 0.9602],
[0.1768, 0.5452, 0.956],
[0.1764, 0.5504, 0.9516],
[0.1755, 0.5554, 0.9473],
[0.174, 0.5605, 0.9432],
[0.1716, 0.5655, 0.9393],
[0.1686, 0.5705, 0.9357],
[0.1649, 0.5755, 0.9323],
[0.161, 0.5805, 0.9289],
[0.1573, 0.5854, 0.9254],
[0.154, 0.5902, 0.9218],
[0.1513, 0.595, 0.9182],
[0.1492, 0.5997, 0.9147],
[0.1475, 0.6043, 0.9113],
[0.1461, 0.6089, 0.908],
[0.1446, 0.6135, 0.905],
[0.1429, 0.618, 0.9022],
[0.1408, 0.6226, 0.8998],
[0.1383, 0.6272, 0.8975],
[0.1354, 0.6317, 0.8953],
[0.1321, 0.6363, 0.8932],
[0.1288, 0.6408, 0.891],
[0.1253, 0.6453, 0.8887],
[0.1219, 0.6497, 0.8862],
[0.1185, 0.6541, 0.8834],
[0.1152, 0.6584, 0.8804],
[0.1119, 0.6627, 0.877],
[0.1085, 0.6669, 0.8734],
[0.1048, 0.671, 0.8695],
[0.1009, 0.675, 0.8653],
[0.0964, 0.6789, 0.8609],
[0.0914, 0.6828, 0.8562],
[0.0855, 0.6865, 0.8513],
[0.0789, 0.6902, 0.8462],
[0.0713, 0.6938, 0.8409],
[0.0628, 0.6972, 0.8355],
[0.0535, 0.7006, 0.8299],
[0.0433, 0.7039, 0.8242],
[0.0328, 0.7071, 0.8183],
[0.0234, 0.7103, 0.8124],
[0.0155, 0.7133, 0.8064],
[0.0091, 0.7163, 0.8003],
[0.0046, 0.7192, 0.7941],
[0.0019, 0.722, 0.7878],
[0.0009, 0.7248, 0.7815],
[0.0018, 0.7275, 0.7752],
[0.0046, 0.7301, 0.7688],
[0.0094, 0.7327, 0.7623],
[0.0162, 0.7352, 0.7558],
[0.0253, 0.7376, 0.7492],
[0.0369, 0.74, 0.7426],
[0.0504, 0.7423, 0.7359],
[0.0638, 0.7446, 0.7292],
[0.077, 0.7468, 0.7224],
[0.0899, 0.7489, 0.7156],
[0.1023, 0.751, 0.7088],
[0.1141, 0.7531, 0.7019],
[0.1252, 0.7552, 0.695],
[0.1354, 0.7572, 0.6881],
[0.1448, 0.7593, 0.6812],
[0.1532, 0.7614, 0.6741],
[0.1609, 0.7635, 0.6671],
[0.1678, 0.7656, 0.6599],
[0.1741, 0.7678, 0.6527],
[0.1799, 0.7699, 0.6454],
[0.1853, 0.7721, 0.6379],
[0.1905, 0.7743, 0.6303],
[0.1954, 0.7765, 0.6225],
[0.2003, 0.7787, 0.6146],
[0.2061, 0.7808, 0.6065],
[0.2118, 0.7828, 0.5983],
[0.2178, 0.7849, 0.5899],
[0.2244, 0.7869, 0.5813],
[0.2318, 0.7887, 0.5725],
[0.2401, 0.7905, 0.5636],
[0.2491, 0.7922, 0.5546],
[0.2589, 0.7937, 0.5454],
[0.2695, 0.7951, 0.536],
[0.2809, 0.7964, 0.5266],
[0.2929, 0.7975, 0.517],
[0.3052, 0.7985, 0.5074],
[0.3176, 0.7994, 0.4975],
[0.3301, 0.8002, 0.4876],
[0.3424, 0.8009, 0.4774],
[0.3548, 0.8016, 0.4669],
[0.3671, 0.8021, 0.4563],
[0.3795, 0.8026, 0.4454],
[0.3921, 0.8029, 0.4344],
[0.405, 0.8031, 0.4233],
[0.4184, 0.803, 0.4122],
[0.4322, 0.8028, 0.4013],
[0.4463, 0.8024, 0.3904],
[0.4608, 0.8018, 0.3797],
[0.4753, 0.8011, 0.3691],
[0.4899, 0.8002, 0.3586],
[0.5044, 0.7993, 0.348],
[0.5187, 0.7982, 0.3374],
[0.5329, 0.797, 0.3267],
[0.547, 0.7957, 0.3159],
[0.5609, 0.7943, 0.305],
[0.5748, 0.7929, 0.2941],
[0.5886, 0.7913, 0.2833],
[0.6024, 0.7896, 0.2726],
[0.6161, 0.7878, 0.2622],
[0.6297, 0.7859, 0.2521],
[0.6433, 0.7839, 0.2423],
[0.6567, 0.7818, 0.2329],
[0.6701, 0.7796, 0.2239],
[0.6833, 0.7773, 0.2155],
[0.6963, 0.775, 0.2075],
[0.7091, 0.7727, 0.1998],
[0.7218, 0.7703, 0.1924],
[0.7344, 0.7679, 0.1852],
[0.7468, 0.7654, 0.1782],
[0.759, 0.7629, 0.1717],
[0.771, 0.7604, 0.1658],
[0.7829, 0.7579, 0.1608],
[0.7945, 0.7554, 0.157],
[0.806, 0.7529, 0.1546],
[0.8172, 0.7505, 0.1535],
[0.8281, 0.7481, 0.1536],
[0.8389, 0.7457, 0.1546],
[0.8495, 0.7435, 0.1564],
[0.86, 0.7413, 0.1587],
[0.8703, 0.7392, 0.1615],
[0.8804, 0.7372, 0.165],
[0.8903, 0.7353, 0.1695],
[0.9, 0.7336, 0.1749],
[0.9093, 0.7321, 0.1815],
[0.9184, 0.7308, 0.189],
[0.9272, 0.7298, 0.1973],
[0.9357, 0.729, 0.2061],
[0.944, 0.7285, 0.2151],
[0.9523, 0.7284, 0.2237],
[0.9606, 0.7285, 0.2312],
[0.9689, 0.7292, 0.2373],
[0.977, 0.7304, 0.2418],
[0.9842, 0.733, 0.2446],
[0.99, 0.7365, 0.2429],
[0.9946, 0.7407, 0.2394],
[0.9966, 0.7458, 0.2351],
[0.9971, 0.7513, 0.2309],
[0.9972, 0.7569, 0.2267],
[0.9971, 0.7626, 0.2224],
[0.9969, 0.7683, 0.2181],
[0.9966, 0.774, 0.2138],
[0.9962, 0.7798, 0.2095],
[0.9957, 0.7856, 0.2053],
[0.9949, 0.7915, 0.2012],
[0.9938, 0.7974, 0.1974],
[0.9923, 0.8034, 0.1939],
[0.9906, 0.8095, 0.1906],
[0.9885, 0.8156, 0.1875],
[0.9861, 0.8218, 0.1846],
[0.9835, 0.828, 0.1817],
[0.9807, 0.8342, 0.1787],
[0.9778, 0.8404, 0.1757],
[0.9748, 0.8467, 0.1726],
[0.972, 0.8529, 0.1695],
[0.9694, 0.8591, 0.1665],
[0.9671, 0.8654, 0.1636],
[0.9651, 0.8716, 0.1608],
[0.9634, 0.8778, 0.1582],
[0.9619, 0.884, 0.1557],
[0.9608, 0.8902, 0.1532],
[0.9601, 0.8963, 0.1507],
[0.9596, 0.9023, 0.148],
[0.9595, 0.9084, 0.145],
[0.9597, 0.9143, 0.1418],
[0.9601, 0.9203, 0.1382],
[0.9608, 0.9262, 0.1344],
[0.9618, 0.932, 0.1304],
[0.9629, 0.9379, 0.1261],
[0.9642, 0.9437, 0.1216],
[0.9657, 0.9494, 0.1168],
[0.9674, 0.9552, 0.1116],
[0.9692, 0.9609, 0.1061],
[0.9711, 0.9667, 0.1001],
[0.973, 0.9724, 0.0938],
[0.9749, 0.9782, 0.0872],
[0.9769, 0.9839, 0.0805]]

parula_map = LinearSegmentedColormap.from_list('parula', cm_data)
matplotlib.rcParams['animation.ffmpeg_path'] = r'C:\ffmpeg\bin\ffmpeg.exe'

def triangle(length, amplitude):
    section = length // 4
    for direction in (1, -1):
        for i in range(section):
            yield i * (amplitude / section) * direction
        for i in range(section):
            yield (amplitude - (i * (amplitude / section))) * direction


def arnold_tongue_save(gamma_0, mu, nu, file):
    nu_positive_grid = np.arange(0, 2, 0.01)
    nu_negative_grid = - np.flip(nu_positive_grid)
    nu_grid = np.append(nu_negative_grid, nu_positive_grid)
    plt.plot(nu_positive_grid, np.sqrt(nu_positive_grid ** 2 + mu ** 2), c='k', linestyle='--')
    plt.fill_between(nu_positive_grid, np.ones(len(nu_positive_grid)) * mu,
                     np.sqrt(nu_positive_grid ** 2 + mu ** 2),
                     facecolor=(92 / 255, 43 / 255, 228 / 255, 0.4))
    plt.plot(nu_negative_grid, np.sqrt(nu_negative_grid ** 2 + mu ** 2), c='k', linestyle='--')
    plt.fill_between(nu_negative_grid, np.ones(len(nu_negative_grid)) * mu,
                     np.sqrt(nu_negative_grid ** 2 + mu ** 2),
                     facecolor=(0, 1, 0, 0.4))
    plt.plot(nu_grid, np.ones(len(nu_grid)) * mu, c='k', linestyle='--')
    plt.fill_between(nu_grid, 2, np.sqrt(nu_grid ** 2 + mu ** 2),
                     facecolor=(1, 0, 0, 0.4))
    plt.fill_between(nu_grid, np.ones(len(nu_grid)) * mu, 0,
                     facecolor=(1, 1, 0, 0.4))
    plt.scatter(nu, gamma_0, c='k', zorder=10)
    plt.title('Arnold Tongue', size='25')
    plt.xlabel('$\\nu$', size='25')
    plt.ylabel('$\gamma$', size='25')
    plt.xlim([-1, 1])
    plt.ylim([0, 1])
    plt.grid(linestyle='--', alpha=0.5)
    plt.savefig(file + '/arnold_tongue.png')
    plt.close()


def arnold_tongue_show(gamma_0, mu, nu):
    nu_positive_grid = np.arange(0, 2, 0.01)
    nu_negative_grid = - np.flip(nu_positive_grid)
    nu_grid = np.append(nu_negative_grid, nu_positive_grid)
    plt.plot(nu_positive_grid, np.sqrt(nu_positive_grid ** 2 + mu ** 2), c='k', linestyle='--')
    plt.fill_between(nu_positive_grid, np.ones(len(nu_positive_grid)) * mu,
                     np.sqrt(nu_positive_grid ** 2 + mu ** 2),
                     facecolor=(92 / 255, 43 / 255, 228 / 255, 0.4))
    plt.plot(nu_negative_grid, np.sqrt(nu_negative_grid ** 2 + mu ** 2), c='k', linestyle='--')
    plt.fill_between(nu_negative_grid, np.ones(len(nu_negative_grid)) * mu,
                     np.sqrt(nu_negative_grid ** 2 + mu ** 2),
                     facecolor=(0, 1, 0, 0.4))
    plt.plot(nu_grid, np.ones(len(nu_grid)) * mu, c='k', linestyle='--')
    plt.fill_between(nu_grid, 2, np.sqrt(nu_grid ** 2 + mu ** 2),
                     facecolor=(1, 0, 0, 0.4))
    plt.fill_between(nu_grid, np.ones(len(nu_grid)) * mu, 0,
                     facecolor=(1, 1, 0, 0.4))
    plt.scatter(nu, gamma_0, c='k', zorder=10)
    plt.title('Arnold Tongue', size='25')
    plt.xlabel('$\\nu$', size='25')
    plt.ylabel('$\gamma$', size='25')
    plt.xlim([-1, 1])
    plt.ylim([0, 1])
    plt.grid(linestyle='--', alpha=0.5)
    plt.show()

def arnold_tongue_exp_show(gamma_0, mu, nu):
    nu_positive_grid = np.arange(0, 2, 0.01)
    nu_negative_grid = - np.flip(nu_positive_grid)
    nu_grid = np.append(nu_negative_grid, nu_positive_grid)
    plt.plot(nu_positive_grid, np.sqrt(nu_positive_grid ** 2 + mu ** 2), c='k', linestyle='--')
    plt.fill_between(nu_positive_grid, np.ones(len(nu_positive_grid)) * mu,
                     np.sqrt(nu_positive_grid ** 2 + mu ** 2),
                     facecolor=(92 / 255, 43 / 255, 228 / 255, 0.4))
    plt.plot(nu_negative_grid, np.sqrt(nu_negative_grid ** 2 + mu ** 2), c='k', linestyle='--')
    plt.fill_between(nu_negative_grid, np.ones(len(nu_negative_grid)) * mu,
                     np.sqrt(nu_negative_grid ** 2 + mu ** 2),
                     facecolor=(0, 1, 0, 0.4))
    plt.plot(nu_grid, np.ones(len(nu_grid)) * mu, c='k', linestyle='--')
    plt.fill_between(nu_grid, 2, np.sqrt(nu_grid ** 2 + mu ** 2),
                     facecolor=(1, 0, 0, 0.4))
    plt.fill_between(nu_grid, np.ones(len(nu_grid)) * mu, 0,
                     facecolor=(1, 1, 0, 0.4))
    plt.scatter(nu, gamma_0, c='k', zorder=10)
    plt.title('Arnold Tongue', size='25')
    plt.xlabel('$\\nu$', size='25')
    plt.ylabel('$\gamma$', size='25')
    plt.xlim([-0.25, 0.25])
    plt.ylim([0, 0.25])
    plt.grid(linestyle='--', alpha=0.5)
    plt.show()


def pix_to_mm(img, scale):
    def click_event(event, x, y, flags, params):
        if event == cv2.EVENT_LBUTTONDOWN:
            points_i = (x, y)
            points.append(points_i)
            cv2.circle(img, (x, y), radius=4, color=(0, 0, 255), thickness=-1)
            cv2.imshow('image', img)
            if len(points) >= 2:
                cv2.line(img, (points[-1]), (points[-2]), (0, 255, 0), thickness=2, lineType=8)
                cv2.circle(img, (points[-1]), radius=4, color=(0, 0, 255), thickness=-1)
                cv2.circle(img, (points[-2]), radius=4, color=(0, 0, 255), thickness=-1)
            cv2.imshow('image', img)

    cv2.imshow('image', img)
    points = []
    cv2.setMouseCallback('image', click_event)
    cv2.waitKey(0)
    pix_to_mm = (40 / np.abs(points[-1][0] - points[-2][0])) * scale
    cv2.destroyAllWindows()
    return pix_to_mm

def injection_length(img, scale):
    def click_event(event, x, y, flags, params):
        if event == cv2.EVENT_LBUTTONDOWN:
            points_i = (x, y)
            points.append(points_i)
            cv2.circle(img, (x, y), radius=4, color=(0, 0, 255), thickness=-1)
            cv2.imshow('image', img)
            if len(points) >= 2:
                cv2.line(img, (points[-1]), (points[-2]), (0, 255, 0), thickness=2, lineType=8)
                cv2.circle(img, (points[-1]), radius=4, color=(0, 0, 255), thickness=-1)
                cv2.circle(img, (points[-2]), radius=4, color=(0, 0, 255), thickness=-1)
            cv2.imshow('image', img)
    cv2.imshow('image', img)
    points = []
    cv2.setMouseCallback('image', click_event)
    cv2.waitKey(0)
    injection_left, injection_right = np.abs(points[-1][0] / scale), np.abs(points[-2][0] / scale)
    cv2.destroyAllWindows()
    return injection_left, injection_right


def define_window(img, scale):
    def click_event(event, x, y, flags, params):
        if event == cv2.EVENT_LBUTTONDOWN:
            points_i = (x, y)
            points.append(points_i)
            cv2.circle(img, (x, y), radius=4, color=(0, 0, 255), thickness=-1)
            cv2.imshow('image', img)
            if len(points) >= 2:
                cv2.line(img, (points[-1]), (points[-2]), (0, 255, 0), thickness=2, lineType=8)
                cv2.circle(img, (points[-1]), radius=4, color=(0, 0, 255), thickness=-1)
                cv2.circle(img, (points[-2]), radius=4, color=(0, 0, 255), thickness=-1)
            cv2.imshow('image', img)

    cv2.imshow('image', img)
    points = []
    cv2.setMouseCallback('image', click_event)
    cv2.waitKey(0)
    window_left_0, window_right_0 = np.abs(points[-1][0] * 1), np.abs(points[-2][0] * 1)
    t_init_pix = np.abs(points[-1][1])
    cv2.destroyAllWindows()
    return window_left_0, window_right_0, t_init_pix

def campos_ligeros(campos, n, Nt, Nx, T):
    t_ligero = np.linspace(0, T, int(Nt / n))
    campos_light = []
    for k in range(len(campos)):
        campo_ligero = np.zeros((int(Nt / n), Nx))
        for i in range(0, len(campos[k][:, 0]) - 1, n):
            campo_ligero[int(i / n), :] = campos[k][i, :]
        campos_light.append(campo_ligero)
    return campos_light, t_ligero


def truncate(num, n):
    integer = int(num * (10**n))/(10**n)
    return float(integer)


def nombre_pndls_LLG(gamma, mu, nu):
    mu_st = str(round(float(mu), 4))
    gamma_st = str(round(float(gamma), 4))
    nu_st = str(round(float(nu), 4))
    nombre = '/gaussian/mu=' + mu_st + '/gamma=' + gamma_st + '/nu=' + nu_st
    return nombre

def nombre_pndls_gaussian(gamma, mu, nu, sigma):
    sigma_st = str(round(float(sigma), 4))
    mu_st = str(round(float(mu), 4))
    gamma_st = str(round(float(gamma), 4))
    nu_st = str(round(float(nu), 4))
    nombre = '/gaussian/mu=' + mu_st + '/gamma=' + gamma_st + '/nu=' + nu_st + '/sigma=' + sigma_st
    return nombre

def nombre_pndls_squared(gamma, mu, nu, length):
    length_st = str(round(float(length), 4))
    mu_st = str(round(float(mu), 4))
    gamma_st = str(round(float(gamma), 4))
    nu_st = str(round(float(nu), 4))
    nombre = '/squared/mu=' + mu_st + '/gamma=' + gamma_st + '/nu=' + nu_st + '/length=' + length_st
    return nombre

def nombre_pndls_bigaussian(gamma, mu, nu, sigma1, sigma2, dist, fase):
    gamma_st = str(truncate(gamma, 4))
    mu_st = str(truncate(mu, 4))
    nu_st = str(truncate(nu, 4))
    sigma1_st = str(truncate(sigma1, 4))
    sigma2_st = str(truncate(sigma2, 4))
    dist_st = str(truncate(dist, 4))
    fase_st = str(truncate(fase / np.pi, 4)) + 'pi'
    nombre = '/bigaussian/mu=' + mu_st + '/gamma=' + gamma_st + '/nu=' + nu_st + '/fase=' + fase_st + '/sigma_1=' + \
             sigma1_st + '_sigma_2=' + sigma2_st + '\\distancia=' + dist_st
    return nombre


def guardar_txt(path, file, **kwargs): # upgradear a diccionario para nombre de variables
    if file == 'no':
        pathfile = path
    else:
        pathfile = path + file
    if os.path.exists(pathfile) == False:
        os.makedirs(pathfile)
    for key, value in kwargs.items():
        np.savetxt(pathfile + '\\' + key + ".txt", value)


def guardar_csv(path, file, **kwargs): # upgradear a diccionario para nombre de variables
    if file == 'no':
        pathfile = path
    else:
        pathfile = path + file
    if os.path.exists(pathfile) == False:
        os.makedirs(pathfile)
    for key, value in kwargs.items():
        np.savetxt(pathfile + '\\' + key + ".csv", value)


def random_transposition(k, N):
    return np.transpose(np.array([k] * N))


def sparse_DD(Nx, dx):
    data = np.ones((3, Nx))
    data[1] = -2 * data[1]
    diags = [-1, 0, 1]
    D2 = sparse.spdiags(data, diags, Nx, Nx) / (dx ** 2)
    D2 = sparse.lil_matrix(D2)
    D2[0, -1] = 1 / (dx ** 2)
    D2[-1, 0] = 1 / (dx ** 2)
    return D2

def sparse_DD_neumann(Nx, dx):
    data = np.ones((3, Nx))
    data[1] = -2 * data[1]
    diags = [-1, 0, 1]
    D2 = sparse.spdiags(data, diags, Nx, Nx) / (dx ** 2)
    D2 = sparse.lil_matrix(D2)
    # Condiciones de borde de Neumann: ajustar la primera y última fila
    D2[0, 0] = -1 / (dx ** 2)
    D2[0, 1] = 1 / (dx ** 2)
    D2[-1, -1] = -1 / (dx ** 2)
    D2[-1, -2] = 1 / (dx ** 2)
    return D2.tocsr()

def sparse_D_periodic(Nx, dx):
    data = np.zeros((3, Nx))
    data[0, :] = -0.5
    data[2, :] = 0.5
    diags = [-1, 0, 1]
    D = sparse.spdiags(data, diags, Nx, Nx, format='lil') / dx
    D[0, -1] = -0.5 / dx
    D[-1, 0] = 0.5 / dx
    return D.tocsr()

def sparse_DD_periodic(Nx, dx):
    data = np.ones((3, Nx))
    data[1, :] = -2.0
    diags = [-1, 0, 1]
    D2 = sparse.spdiags(data, diags, Nx, Nx, format='lil') / dx**2
    D2[0, -1] = 1.0 / dx**2
    D2[-1, 0] = 1.0 / dx**2
    return D2.tocsr()

def sparse_DD_absorbing(Nx, dx, absorption_coeff):
    data = np.ones((3, Nx))
    data[1] = -2 * data[1]
    diags = [-1, 0, 1]
    D2 = sparse.spdiags(data, diags, Nx, Nx) / (dx ** 2)
    D2 = sparse.lil_matrix(D2)

    # Apply absorbing boundary conditions using exponential damping
    sigma_L = absorption_coeff / (dx ** 2)
    sigma_R = absorption_coeff / (dx ** 2)

    # Left boundary
    D2[0, 0] = -2 / (dx ** 2) - sigma_L
    D2[0, 1] = 2 / (dx ** 2)

    # Right boundary
    D2[-1, -1] = -2 / (dx ** 2) - sigma_R
    D2[-1, -2] = 2 / (dx ** 2)

    return D2.tocsr()

def sparse_D_neumann(Nx, dx):
    data = np.zeros((3, Nx))
    data[0, 1:] = -0.5  # lower
    data[2, :-1] = 0.5  # upper
    diags = [-1, 0, 1]
    D = sparse.spdiags(data, diags, Nx, Nx) / dx
    D = D.tolil()

    # Neumann BCs: forward/backward difference at the ends
    D[0, 0] = -1 / dx
    D[0, 1] = 1 / dx
    D[-1, -2] = -1 / dx
    D[-1, -1] = 1 / dx
    return D.tocsr()

    return D.tocsr()

def sparse_D(Nx, dx):
    data = np.ones((3, Nx))
    data[0] = - 1 * data[0]
    data[1] = 0 * data[1]
    data[2] = 1 * data[2]
    diags = [-1, 0, 1]
    D1 = sparse.spdiags(data, diags, Nx, Nx) / (2 * dx)
    D1 = sparse.lil_matrix(D1)
    D1[0, 0] = - 1.0 / (2 * dx)
    D1[-1, -1] = 0.0 / (2 * dx)
    return D1

def sparse_D_neumann_4order(Nx, dx):
    data = np.ones((5, Nx))
    data[0] = 1 * data[0]
    data[1] = - 8 * data[1]
    data[2] = 0 * data[2]
    data[3] = 8 * data[3]
    data[4] = - 1 * data[4]
    diags = [-2, -1, 0, 1, 2]
    D1 = sparse.spdiags(data, diags, Nx, Nx) / (12 * dx)
    D1 = sparse.lil_matrix(D1)
    D1[0, 0] = 0
    D1[-1, -1] = 0
    D1[1, 0] = 0
    D1[-2, -1] = 0
    return D1

def sparse_DDD(Nx, dx):
    data = np.ones((5, Nx))
    data[0] = - 0.5 * data[0]
    data[1] = data[1]
    data[2] = 0 * data[2]
    data[3] = - 1 * data[3]
    data[4] = 0.5 * data[4]
    diags = [-2, -1, 0, 1, 2]
    D3 = sparse.spdiags(data, diags, Nx, Nx) / (dx ** 3)
    D3 = sparse.lil_matrix(D3)
    D3[-1, 0] = - 1 / (dx ** 3)
    D3[-1, 1] = 0.5 / (dx ** 3)
    D3[-2, 0] = 0.5 / (dx ** 3)

    D3[0, -1] = 1 / (dx ** 3)
    D3[0, -2] = - 0.5 / (dx ** 3)
    D3[1, -1] = - 0.5 / (dx ** 3)
    return D3


def  sparse_DDDD_neumann(Nx, dx):
    # Interior stencil for 4th derivative: [1, -4, 6, -4, 1]
    diagonals = [
        np.ones(Nx - 2),                # -2
        -4 * np.ones(Nx - 1),           # -1
        6 * np.ones(Nx),                #  0
        -4 * np.ones(Nx - 1),           # +1
        np.ones(Nx - 2)                 # +2
    ]
    offsets = [-2, -1, 0, 1, 2]

    D4 = sparse.diags(diagonals, offsets, shape=(Nx, Nx), format='lil') / dx**4

    # Neumann boundary conditions (approximate with forward/backward biased stencils)
    # Left boundary (x = 0)
    D4[0, 0:5] = [16, -64, 96, -64, 16]  # From 5-point forward difference
    D4[1, 0:5] = [1, -4, 6, -4, 1]       # Normal central stencil starts at i=1

    # Right boundary (x = Nx-1)
    D4[-1, -5:] = [16, -64, 96, -64, 16]  # From 5-point backward difference
    D4[-2, -5:] = [1, -4, 6, -4, 1]       # Normal central stencil at i=Nx-2

    return D4.tocsr()


def sparse_DDDD(Nx, dx):
    data = np.ones((5, Nx))
    data[1] = -4 * data[1]
    data[4] = 6 * data[2]
    data[3] = -4 * data[3]
    diags = [-2, -1, 0, 1, 2]
    D4 = sparse.spdiags(data, diags, Nx, Nx) / (dx ** 4)
    D4 = sparse.lil_matrix(D4)
    D4[0, -1] = -4 / (dx ** 4)
    D4[0, -2] = 1 / (dx ** 4)
    D4[1, -1] = 1 / (dx ** 4)

    D4[-1, 0] = -4 / (dx ** 4)
    D4[-1, 1] = 1 / (dx ** 4)
    D4[-2, 0] = 1 / (dx ** 4)
    return D4


def sparse_id(Nx):
    data = np.ones((1, Nx))
    diags = [0]
    id = sparse.spdiags(data, diags, Nx, Nx)
    id = sparse.lil_matrix(id)
    return id


def Dxx(DD, f):
    dd_f = DD @ f
    return dd_f

def Dx(D, f):
    d_f = D @ f
    return d_f


def Der(D, f):
    d_f = D @ f
    return d_f

def Dxxx(DDD, f):
    ddd_f = DDD @ f
    return ddd_f


def filtro_array(n, funcion):
    # the larger n is, the smoother curve will be
    b = [1.0 / n] * n
    a = 1
    Z_filtered = filtfilt(b, a, funcion)
    return Z_filtered

#COMENTARIO

def filtro_superficie(Z, intensidad, sentido):
    X_len = len(Z[:, 0])
    Y_len = len(Z[0, :])
    FILT = np.zeros((X_len, Y_len))
    if sentido == 'X':
        for i in range(X_len):
            filtered = filtro_array(intensidad, Z[i, :])
            FILT[i, :] = filtered
    elif sentido == 'Y':
        for i in range(Y_len):
            filtered = filtro_array(intensidad, Z[:, i])
            FILT[:, i] = filtered
    elif sentido == 'XY':
        for i in range(X_len):
            filtered = filtro_array(intensidad, Z[i, :])
            FILT[i, :] = filtered
        for i in range(Y_len):
            filtered = filtro_array(intensidad, FILT[:, i])
            FILT[:, i] = filtered
    elif sentido == 'YX':
        for i in range(Y_len):
            filtered = filtro_array(intensidad, Z[:, i])
            FILT[:, i] = filtered
        for i in range(X_len):
            filtered = filtro_array(intensidad, FILT[i, :])
            FILT[i, :] = filtered
    return FILT

def pdnlS_name(parameters, order):
    rounded_params = [np.round(param, 3) for param in parameters]
    str_params = [str(param) for param in rounded_params]
    str_params_final = []
    for str_param in str_params:
        splited_params = str_param.split(".")
        len_digits = len(splited_params[-1])
        if len_digits == 1 and len(splited_params) != 1:
            str_param = str_param + "00"
        elif len_digits == 2 and len(splited_params) != 1:
            str_param = str_param + "0"
        elif len(splited_params) == 1:
            str_param = str_param + ".000"
        str_params_final.append(str_param)
    if order == "ABMNSG":
        directory_name = "/alpha=" + str_params_final[0] + "/beta=" + str_params_final[1] + "/mu=" + str_params_final[2] + "/nu=" + str_params_final[3] + "/sigma=" + str_params_final[4] + "/gamma=" + str_params_final[5]
    elif order == "ABMNGS":
        directory_name = "/alpha=" + str_params_final[0] + "/beta=" + str_params_final[1] + "/mu=" + str_params_final[2] + "/nu=" + str_params_final[3] + "/gamma=" + str_params_final[5] + "/sigma=" + str_params_final[4]
    return directory_name

def pdnlS_bigauss_name(parameters, order):
    rounded_params = [np.round(param, 4) for param in parameters]
    str_params = [str(param) for param in rounded_params]
    str_params_final = []
    for str_param in str_params:
        splited_params = str_param.split(".")
        len_digits = len(splited_params[-1])
        if len_digits == 1 and len(splited_params) != 1:
            str_param = str_param + "000"
        elif len_digits == 2 and len(splited_params) != 1:
            str_param = str_param + "00"
        elif len_digits == 3 and len(splited_params) != 1:
            str_param = str_param + "0"
        elif len(splited_params) == 1:
            str_param = str_param + ".000"
        str_params_final.append(str_param)
    if order == "ABMNSGD":
        directory_name = "/alpha=" + str_params_final[0] + "/beta=" + str_params_final[1] + "/mu=" + str_params_final[2] + "/nu=" + str_params_final[3] + "/sigma=" + str_params_final[4] + "/gamma=" + str_params_final[5] + "/dist=" + str_params_final[6]
    elif order == "ABMNGSD":
        directory_name = "/alpha=" + str_params_final[0] + "/beta=" + str_params_final[1] + "/mu=" + str_params_final[2] + "/nu=" + str_params_final[3] + "/gamma=" + str_params_final[5] + "/sigma=" + str_params_final[4] + "/dist=" + str_params_final[6]
    elif order == "exp_sim":
        directory_name = "/n_pist=" + str_params_final[0] + "/freq=" + str_params_final[1] + "/ampd=" + str_params_final[2] + "/d_pist=" + str_params_final[3]
    return directory_name


def pdnlS_str_parameters(parameters, order):
    rounded_params = [np.round(param, 3) for param in parameters]
    str_params = [str(param) for param in rounded_params]
    str_params_final = []
    for str_param in str_params:
        splited_params = str_param.split(".")
        len_digits = len(splited_params[-1])
        if len_digits == 1 and len(splited_params) != 1:
            str_param = str_param + "00"
        elif len_digits == 2 and len(splited_params) != 1:
            str_param = str_param + "0"
        elif len(splited_params) == 1:
            str_param = str_param + ".000"
        str_params_final.append(str_param)
    return str_params_final



def jacobian_func(U_1, U_2, parameters):
    [alpha, beta, gamma, mu, nu] = parameters
    gamma_1 = gamma[0]
    gamma_2 = gamma[1]

    j_11 = [gamma_1 - mu + 2 * beta * U_1 * U_2]
    j_12 = [nu + gamma_2 + beta * (U_1 ** 2 + U_2 ** 2) + 2 * beta * U_2 ** 2]
    j_21 = [nu + gamma_2 + beta * (U_1 ** 2 + U_2 ** 2) + 2 * beta * U_2 ** 2]
    j_22 = [- nu + gamma_2 - beta * (U_1 ** 2 + U_2 ** 2) - 2 * beta * U_1 ** 2]
    return [j_11, j_12, j_21, j_22]


def jacobian_pdnlS(fields, parameters, operators):
    U_1 = fields[0]
    U_2 = fields[1]
    [alpha, beta, gamma, mu, nu] = parameters
    gamma_1 = gamma[0]
    DD = operators[0]

    j_11 = [gamma_1 - mu + 2 * beta * U_1 * U_2]
    j_12 = [nu + beta * U_1 ** 2 + 3 * beta * U_2 ** 2]
    j_21 = [- nu - beta * U_1 ** 2 - 3 * beta * U_1 ** 2]
    j_22 = [- (gamma_1 + mu) - 2 * beta * U_2 * U_1]

    diagonals = j_11
    J_11 = diags(diagonals, [0])
    J_11 = J_11.toarray()

    diagonals = j_12
    J_12 = diags(diagonals, [0]) + alpha * DD
    J_12 = J_12.toarray()

    diagonals = j_21
    J_21 = diags(diagonals, [0]) - alpha * DD
    J_21 = J_21.toarray()

    diagonals = j_22
    J_22 = diags(diagonals, [0])
    J_22 = J_22.toarray()

    J_1 = np.concatenate((J_11, J_12), axis=1)
    J_2 = np.concatenate((J_21, J_22), axis=1)
    J = np.concatenate((J_1, J_2), axis=0)
    return J

def jacobian_pdnlS_spectral(fields, parameters, kappa, Nx):
    U_1 = fields[0]
    U_2 = fields[1]
    Uhat_1 = np.fft.fft(U_1)
    Uhat_2 = np.fft.fft(U_2)
    dd_Uhat_1 = -np.power(kappa, 2) * Uhat_1
    dd_Uhat_2 = -np.power(kappa, 2) * Uhat_2
    [alpha, beta, gamma, mu, nu] = parameters
    gamma_1 = gamma[0]

    J_11_a = gamma_1 - mu + 2 * beta * U_1 * U_2
    J_12_a = nu + beta * U_1 ** 2 + 3 * beta * U_2 ** 2
    J_21_a = - nu - beta * U_1 ** 2 - 3 * beta * U_1 ** 2
    J_22_a = - (gamma_1 + mu) - 2 * beta * U_2 * U_1

    J_11_b = np.zeros(Nx)
    J_12_b = alpha
    J_21_b = - alpha
    J_22_b = np.zeros(Nx)

    J_1_a = np.concatenate((J_11_a, J_12_a), axis=1)
    J_2_a = np.concatenate((J_21_a, J_22_a), axis=1)
    J_a = np.concatenate((J_1_a, J_2_a), axis=0)

    J_1_b = np.concatenate((J_11_b, J_12_b), axis=1)
    J_2_b = np.concatenate((J_21_b, J_22_b), axis=1)
    J_b = np.concatenate((J_1_a, J_2_a), axis=0)
    return np.array([J_a, J_b])


def jacobian_wave(fields, parameters, operators, Nx):
    U_1 = fields[0]
    U_2 = fields[1]
    [c] = parameters
    DD = operators[0]

    j_11 = [np.zeros(Nx)]
    j_12 = [np.ones(Nx)]
    j_21 = [np.zeros(Nx)]
    j_22 = [np.zeros(Nx)]

    diagonals = j_11
    J_11 = diags(diagonals, [0])
    J_11 = J_11.toarray()

    diagonals = j_12
    J_12 = diags(diagonals, [0])
    J_12 = J_12.toarray()

    diagonals = j_21
    J_21 = diags(diagonals, [0]) + c ** 2 * DD
    J_21 = J_21.toarray()

    diagonals = j_22
    J_22 = diags(diagonals, [0])
    J_22 = J_22.toarray()

    J_1 = np.concatenate((J_11, J_12), axis=0)
    J_2 = np.concatenate((J_21, J_22), axis=0)
    J = np.concatenate((J_1, J_2), axis=1)
    return J

def jacobian_lug_lef(fields, parameters, operators, Nx):
    U_1 = fields[0]
    U_2 = fields[1]
    delta = parameters[0][0]
    DD = operators[0]

    j_11 = [-1 - 2 * U_1 * U_2]
    j_12 = [delta - U_1 ** 2 - 3 * U_2 ** 2]
    j_21 = [-delta + U_2 ** 2 + 3 * U_1 ** 2]
    j_22 = [-1 + 2 * U_2 * U_1]

    diagonals = j_11
    J_11 = diags(diagonals, [0])
    J_11 = J_11.toarray()
    #J_11[0, 0] = 100

    diagonals = j_12
    J_12 = diags(diagonals, [0]) - DD
    J_12 = J_12.toarray()
    #J_12[0, 0] = 200

    diagonals = j_21
    J_21 = diags(diagonals, [0]) + DD
    J_21 = J_21.toarray()
    #J_21[0, 0] = 300

    diagonals = j_22
    J_22 = diags(diagonals, [0])
    J_22 = J_22.toarray()
    #J_22[0, 0] = 400

    J_1 = np.concatenate((J_11, J_12), axis=1)
    J_2 = np.concatenate((J_21, J_22), axis=1)
    J = np.concatenate((J_1, J_2), axis=0)
    return J


def jacobian_KS(fields, parameters, operators, Nx):
    U_1 = fields[0]
    D1 = operators[0]
    D2 = operators[1]
    D4 = operators[2]

    j = [- 2 * D1 @ U_1]

    diagonals = j
    J = diags(diagonals, [0]) @ D1 - D4 - D2
    J = J.toarray()
    return J

def time_propagator(type, I, J, Q, dt):
    if type == "I_Jdt_RK4":
        k_1 = J @ Q  # np.exp(J * dt) @ Q #
        k_2 = J @ (Q + 0.5 * dt * k_1)  # np.exp(J * dt) @ (Q + 0.5 * dt * k_1) #
        k_3 = J @ (Q + 0.5 * dt * k_2)  # np.exp(J * dt) @ (Q + 0.5 * dt * k_2) #
        k_4 = J @ (Q + dt * k_3)  # np.exp(J * dt) @ (Q + dt * k_3) #
        U_new = Q + dt * (k_1 + 2 * k_2 + 2 * k_3 + k_4) / 6
    elif type == "exp":
        U_new = np.exp(J * dt) @ Q
    elif type == "I_Jdt":
        U_new = (I + J * dt) @ Q
    return U_new


def fluid_pdnls_parameters(f_i, a_ang, d):
    g = 9790
    l_y = 14.2
    w = 2 * np.pi * (f_i / 2)
    k_y = np.pi / l_y
    k = k_y
    tau = np.tanh(k * d)
    w_1 = np.sqrt(g * k * tau)
    f_0 = w_1 / np.pi
    gamma = f_i ** 2 * a_ang * 8.401093 * 10 ** (-5)
    alpha = (1 / (4 * k ** 2)) * (1 + k * d * ((1 - tau ** 2) / tau))
    #alpha = (1 / (4 * tau)) * (1 - tau ** 2 + tau/k)# término difusivo
    beta = (k ** 2 / 64) * (6 * tau ** 2 - 5 + 16 * tau ** (-2) - 9 * tau ** (-4))  # término no lineal
    nu = 0.5 * ((w / w_1) ** 2 - 1)
    return alpha, beta, nu, gamma, f_0


def max_finder(Z, t_grid, Nt, dt):
    D = sparse_D_neumann_4order(Nt, dt)
    DD = sparse_DD_neumann(Nt, dt)

    #print(np.transpose(Z).shape)
    #print(Nt)

    D1_Z = D.dot(np.transpose(Z))
    D2_Z = DD.dot(np.transpose(Z))

    tau_L_points_max = []
    Z_points_max = []
    I_max = []
    #for i in range(Nt):
    #    if np.sign(D1_Z[i]) != np.sign(D1_Z[i - 1]) and D2_Z[i] < 0 and i != 0 and np.sign(D1_Z[i]) + np.sign(
    #            D1_Z[i - 1]) != 1:
    #        tau_L_points_max.append(t_grid[i])
    #        Z_points_max.append(Z[i])
    #        I_max.append(i)
    Z_points_max = np.array(Z[signal.argrelextrema(Z, np.greater)[0]])
    tau_L_points_max = np.array(t_grid[signal.argrelextrema(Z, np.greater)[0]])
    I_max = np.array(signal.argrelextrema(Z, np.greater)[0])
    return Z_points_max, tau_L_points_max, I_max

def min_finder(Z, t_grid, Nt, dt):
    D = sparse_D_neumann_4order(Nt, dt)
    DD = sparse_DD_neumann(Nt, dt)
    #Z = filtro_array(10, Z)
    #print(np.transpose(Z).shape)
    #print(Nt)

    D1_Z = D.dot(np.transpose(Z))
    D2_Z = DD.dot(np.transpose(Z))

    tau_L_points_max = []
    Z_points_max = []
    I_max = []
    for i in range(Nt):
        if np.sign(D1_Z[i]) != np.sign(D1_Z[i - 1]) and D2_Z[i] > 0 and i != 0 and np.sign(D1_Z[i]) + np.sign(
                D1_Z[i - 1]) != 1:
            tau_L_points_max.append(t_grid[i])
            Z_points_max.append(Z[i])
            I_max.append(i)
    Z_points_max = np.array(Z_points_max)
    tau_L_points_max = np.array(tau_L_points_max)
    I_max = np.array(I_max)
    return Z_points_max, tau_L_points_max, I_max

def phis(alpha, beta, nu, mu, gamma, sigma, X, Y, x_grid, dx):
    gamma_01 = gamma * np.exp(- X ** 2 / (2 * sigma ** 2))
    gamma_02 = gamma * np.exp(- Y ** 2 / (2 * sigma ** 2))
    delta_01 = - nu + np.sqrt(gamma_01 ** 2 - mu ** 2)
    delta_02 = - nu + np.sqrt(gamma_02 ** 2 - mu ** 2)
    theta_01 = 0.5 * np.arccos(mu / gamma_01)
    theta_02 = 0.5 * np.arccos(mu / gamma_02)
    #phi_01 = (1 / (np.cosh(np.sqrt(delta_01 / alpha) * (x_grid - X)))) * np.sqrt(2 * delta_01) * np.cos(theta_01) * (beta ** (-2))
    #phi_02 = - (1 / (np.cosh(np.sqrt(delta_02 / alpha) * (x_grid - Y)))) * np.sqrt(2 * delta_02) * np.sin(theta_02) * (beta ** (-2))
    phi_01 = np.exp(-0.89 * np.sqrt(delta_01 / alpha) ** 2 * (x_grid - X) ** 2) * np.sqrt(2 * delta_01) * np.cos(theta_01) * (beta ** (-2))
    phi_02 = - np.exp(-0.89*np.sqrt(delta_02 / alpha) ** 2 * (x_grid - Y) ** 2) * np.sqrt(2 * delta_02) * np.sin(theta_02) * (beta ** (-2))
    Dphi_01 = np.append(np.diff(phi_01) / dx, 0)
    Dphi_02 = np.append(np.diff(phi_02) / dx, 0)
    return [phi_01, phi_02, Dphi_01, Dphi_02]