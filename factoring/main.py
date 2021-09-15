import numpy as np
from scipy import constants
import matplotlib.pyplot as plt
from number_of_bits import number_of_bits_for_x, number_of_bits_for_y

I = np.eye(2)
X = np.array([[0, 1],
              [1, 0]])
Z = np.array([[1, 0],
              [0,-1]])

IN = (I - X) / 2

def tensor_product(M1, *args):
    if len(args) == 0:
        return M1
    return np.kron(M1, tensor_product(*args))

def apply_to_bit(M, i, n):
    pre = [I] * i
    pos = [I] * ((n - 1) - i)
    sequence = pre + [M] + pos
    return tensor_product(*sequence)

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
        partial_x += (2**(l+1) * (I_nx - apply_to_bit(Z, l, nx)) / 2)

    partial_y = I_ny
    for m in range(ny):
        partial_y += (2**(m+1) * (I_ny - apply_to_bit(Z, m, ny)) / 2)

    HP = N * I_n
    HP = HP - tensor_product(partial_x, partial_y)
    HP = HP ** 2

    HB = np.zeros((2**n, 2**n))
    for i in range(n):
        HB += apply_to_bit((I - X) / 2, i, n)

    B = lambda s: np.sin(np.pi/2 * (np.sin(np.pi * s / 2) ** 2)) ** 2
    A = lambda s: 1 - B(s)
    H = lambda s: (A(s) * HB) + (B(s) * HP)

    t = 0
    ts = []
    eigenvalues = []
    while t <= 1:
        eigenvalues.append(np.linalg.eigvalsh(H(t)))
        ts.append(t)
        t += 0.01

    eigenvalues = np.stack(eigenvalues, axis=1)
    for i, eigenvalue in enumerate(eigenvalues):
        plt.plot(ts, eigenvalue, c="r" if i == 0 else "b")
    plt.xlabel("Tempo")
    plt.ylabel("Energia")
    plt.show()

if __name__ == "__main__":
    main()