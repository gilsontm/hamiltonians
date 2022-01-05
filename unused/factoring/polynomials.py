import ipdb
import __init__
import numpy as np
from number_of_bits import number_of_bits_for_x, number_of_bits_for_y
import math

def set_product(poly1, poly2):
    terms = set()
    for term1 in poly1:
        for term2 in poly2:
            terms.add(term1.union(term2))
    return terms

def variables(N):
    nx = int(number_of_bits_for_x(N))
    ny = int(number_of_bits_for_y(N))
    x = set([frozenset([i])    for i in range(nx)] + [frozenset()])
    y = set([frozenset([i+nx]) for i in range(ny)] + [frozenset()])
    terms = set_product(x, y)
    product = set_product(terms, terms)
    return len(list(filter(lambda s: len(s) >= 3, product)))

def n_bits_for_x(N):
    floor = math.floor(math.sqrt(N))
    floor -= 1 - (floor % 2)
    nx = math.ceil(math.log2(floor)) - 1
    return nx

def n_bits_for_y(N):
    ny = math.ceil(math.log2(math.floor(N / 3))) - 1
    return ny

def main():
    x = np.array(range(2, 50))
    y = [variables(2**n) for n in x]

    n = list(range(2, 1024))
    nx = [n_bits_for_x(2**i) for i in n]
    ny = [n_bits_for_y(2**i) for i in n]

    s = [(nxi * nyi) * (nxi * nyi + nxi + nyi - 3) / 4 for nxi, nyi in zip(nx, ny)]

    # poly3 = x ** 3
    # poly4 = x ** 4
    import matplotlib.pyplot as plt
    plt.plot(n, s, "b", label="Variáveis extras (algébrico)")
    plt.plot(x, y, "r", label="Variáveis extras (numérico)")
    # plt.plot(x, poly3, "g^-", label="n^3")
    # plt.plot(x, poly4, "bs-", label="n^4")
    plt.legend()
    plt.xlabel("Bits de N")
    # plt.ticklabel_format(style='plain')
    plt.show()

if __name__ == "__main__":
    main()