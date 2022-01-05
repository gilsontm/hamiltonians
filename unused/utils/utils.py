import numpy as np
import matplotlib.pyplot as plt


I = np.eye(2)
X = np.array([[0, 1],
              [1, 0]])
Z = np.array([[1, 0],
              [0,-1]])

def tensor_product(M, *args):
    if len(args) == 0:
        return M
    return np.kron(M, tensor_product(*args))

def apply_to_bit(M, i, n):
    pre = [I] * i
    pos = [I] * ((n - 1) - i)
    sequence = pre + [M] + pos
    return tensor_product(*sequence)

def evolve(HB, HP, f = lambda s: s):
    H = lambda s: ((1 - f(s)) * HB) + (f(s) * HP)

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
