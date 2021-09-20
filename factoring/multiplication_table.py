import __init__
from utils.utils import *


M = (I - Z) / 2

def complete():
    """
        |p1p2q1q2|

        (p1 + q1 - 1)**2 + (p2 + q2 - 1)**2 + (p2 * q1 + p1 * q2 - 1)**2 = 0
    """
    H1  =      tensor_product(M, I, I, I)
    H1 +=      tensor_product(I, I, M, I)
    H1 += -1 * tensor_product(I, I, I, I)
    H1  = H1 ** 2

    H2  =      tensor_product(I, M, I, I)
    H2 +=      tensor_product(I, I, I, M)
    H2 += -1 * tensor_product(I, I, I, I)
    H2  = H2 ** 2

    H3  =      tensor_product(I, M, M, I)
    H3 +=      tensor_product(M, I, I, M)
    H3 += -1 * tensor_product(I, I, I, I)
    H3  = H3 ** 2

    HP  = H1 + H2 + H3
    f   = lambda p1, p2, q1, q2: (p1 + q1 - 1)**2 + (p2 + q2 - 1)**2 + (p2 * q1 + p1 * q2 - 1)**2
    return HP, f


def simplified():
    """
        |p1p2q1q2|

        3 - p2 * q1 - p1 * q2 + 2 * p2 * q2 - p2 - q2 + 2 * p1 * q1 - p1 - q1 = 0
    """
    HP  =  3 * tensor_product(I, I, I, I)
    HP += -1 * tensor_product(I, M, M, I)
    HP += -1 * tensor_product(M, I, I, M)
    HP +=  2 * tensor_product(I, M, I, M)
    HP += -1 * tensor_product(I, M, I, I)
    HP += -1 * tensor_product(I, I, I, M)
    HP +=  2 * tensor_product(M, I, M, I)
    HP += -1 * tensor_product(M, I, I, I)
    HP += -1 * tensor_product(I, I, M, I)
    f   = lambda p1, p2, q1, q2: 3 - p2 * q1 - p1 * q2 + 2 * p2 * q2 - p2 - q2 + 2 * p1 * q1 - p1 - q1
    # f   = lambda p1, p2, q1, q2: 3 - p2 * q1 - p1 * q2 - p2 - q2 - p1 - q1
    return HP, f

def more_simplified():
    """
        |p1p2q1q2|

        3 - p2 * q1 - p1 * q2 - p2 - q2 - p1 - q1 = 0
    """
    HP  =  3 * tensor_product(I, I, I, I)
    HP += -1 * tensor_product(I, M, M, I)
    HP += -1 * tensor_product(M, I, I, M)
    # HP +=  2 * tensor_product(I, M, I, M)
    HP += -1 * tensor_product(I, M, I, I)
    HP += -1 * tensor_product(I, I, I, M)
    # HP +=  2 * tensor_product(M, I, M, I)
    HP += -1 * tensor_product(M, I, I, I)
    HP += -1 * tensor_product(I, I, M, I)
    f   = lambda p1, p2, q1, q2: 3 - p2 * q1 - p1 * q2 - p2 - q2 - p1 - q1
    return HP, f

def main():
    H1, f1 = complete()
    H2, f2 = simplified()
    H3, f3 = more_simplified()

    n = 4
    r1 = []
    r2 = []
    r3 = []
    values = []
    for i in range(2**n):
        binary = f"{i:0{n}b}"
        values.append(binary)
        binary = list(map(int, binary))
        r1.append(f1(*binary))
        r2.append(f2(*binary))
        r3.append(f3(*binary))
        print(f"{i:0{n}b}", f1(*binary), f2(*binary), f3(*binary))

    print(H1); print(H2);  print(H3)

if __name__ == "__main__":
    main()
