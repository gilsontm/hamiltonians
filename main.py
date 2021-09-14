import numpy as np
from ket import *
import matplotlib.pyplot as plt
from numpy.lib import NumpyVersion


I = np.eye(2)
X = np.array([[0, 1],
              [1, 0]])
Z = np.array([[1, 0],
              [0,-1]])

P1 =      (I - Z) / 2        # punishes ket1
P0 = (I - (I - Z) / 2)       # punishes ket0

IN = (I - X) / 2

def tensor_product(M1, *args):
    if len(args) == 0:
        return M1
    return np.kron(M1, tensor_product(*args))

def plot(H, t):
    e = np.linalg.eigvalsh(H)
    t = len(e) * [t]
    plt.scatter(t, e, s=1, c="r")

def main():
                                     # 0 1 2 3
                                     # -------
    ' (¬q0 v ¬q1 v        q3) &&  '  # 1 1 _ 0        = 11_0
    ' (¬q0 v       ¬q2 v ¬q3) &&  '  # 1 _ 1 1        = 1_11
    ' (¬q0 v  q1 v  q2      ) &&  '  # 1 0 0 _        = 100_
    ' (¬q0 v  q1 v        q3) &&  '  # 1 0 _ 0        = 10_0
    ' ( q0 v ¬q1 v       ¬q3) &&  '  # 0 1 _ 1        = 01_1
    ' ( q0 v ¬q1 v        q3) &&  '  # 0 1 _ 0        = 01_0
    ' ( q0 v  q1 v ¬q2      ) &&  '  # 0 0 1 _        = 001_
    ' ( q0 v  q1 v  q2      )     '  # 0 0 0 _        = 000_

    I = np.eye(2)

    #                    q0  q1  q2  q3
    HB = (tensor_product(IN,  I,  I,  I) +
           tensor_product( I, IN,  I,  I) +
           tensor_product( I,  I, IN,  I) +
           tensor_product( I,  I,  I, IN))

    #                    q0  q1  q2  q3
    HP = (tensor_product(P1, P1,  I, P0) +
          tensor_product(P1,  I, P1, P1) +
          tensor_product(P1, P0, P0,  I) +
          tensor_product(P1, P0,  I, P0) +
          tensor_product(P0, P1,  I, P1) +
          tensor_product(P0, P1,  I, P0) +
          tensor_product(P0, P0, P1,  I) +
          tensor_product(P0, P0, P0,  I))

    H = lambda s: ((1-s) * HB) + (s * HP)

    t = 0
    while t <= 1:
        plot(H(t), t)
        t += 0.01
    plt.show()

if __name__ == "__main__":
    main()
