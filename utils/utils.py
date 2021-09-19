import numpy as np

I = np.eye(2)
X = np.array([[0, 1],
              [1, 0]])
Z = np.array([[1, 0],
              [0,-1]])

IN = (I - X) / 2

def tensor_product(M, *args):
    if len(args) == 0:
        return M
    return np.kron(M, tensor_product(*args))

def apply_to_bit(M, i, n):
    pre = [I] * i
    pos = [I] * ((n - 1) - i)
    sequence = pre + [M] + pos
    return tensor_product(*sequence)
