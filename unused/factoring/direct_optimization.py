import __init__
import numpy as np
from utils.utils import *
import matplotlib.pyplot as plt
from number_of_bits import number_of_bits_for_x, number_of_bits_for_y


def main():
    N = 21
    nx = int(number_of_bits_for_x(N))
    ny = int(number_of_bits_for_y(N))
    nx = int(nx)
    ny = int(ny)
    n  = nx + ny

    I_nx = np.eye(2**nx)
    I_ny = np.eye(2**ny)
    I_n  = np.eye(2**n)

    partial_x = I_nx
    for l in range(nx):
        partial_x += (2**(l+1) * apply_to_bit((I - Z) / 2, l, nx))

    partial_y = I_ny
    for m in range(ny):
        partial_y += (2**(m+1) * apply_to_bit((I - Z) / 2, m, ny))

    HP = N * I_n
    HP = HP - tensor_product(partial_x, partial_y)
    HP = HP ** 2

    HB = np.zeros((2**n, 2**n))
    for i in range(n):
        HB += apply_to_bit((I - X) / 2, i, n)

    f = lambda s: np.sin(np.pi/2 * (np.sin(np.pi * s / 2) ** 2)) ** 2

    evolve(HB, HP, f)

if __name__ == "__main__":
    main()
