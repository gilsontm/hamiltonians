import numpy as np

def number_of_bits_for_x(N):
    floor = np.floor(np.sqrt(N))
    floor -= 1 - (floor % 2)
    nx = np.ceil(np.log2(floor)) - 1
    return nx

def number_of_bits_for_y(N):
    ny = np.ceil(np.log2(np.floor(N / 3))) - 1
    return ny
