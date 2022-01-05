import __init__
import numpy as np
from utils.utils import *

P1 =      (I - Z) / 2        # punishes ket1
P0 = (I - (I - Z) / 2)       # punishes ket0

IN = (I - X) / 2

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

    evolve(HB, HP)


def other():

    HP1 = (tensor_product(P0, P0, P0,  I) +
           tensor_product(P1,  I, P0, P0) +
           tensor_product( I, P0, P1, P1))

    HP2 = (  tensor_product( I,  I,  I,  I) # part 1
           - tensor_product(P1,  I,  I,  I)
           - tensor_product( I, P1,  I,  I)
           + tensor_product(P1, P1,  I,  I)
           - tensor_product( I,  I, P1,  I)
           + tensor_product(P1,  I, P1,  I)
           + tensor_product( I, P1, P1,  I)
           - tensor_product(P1, P1, P1,  I)

           + tensor_product(P1,  I,  I,  I) # part 2
           - tensor_product(P1,  I, P1,  I)
           - tensor_product(P1,  I,  I, P1)
           + tensor_product(P1,  I, P1, P1)

           + tensor_product( I,  I, P1, P1) # part 3
           - tensor_product( I, P1, P1, P1))

    import ipdb; ipdb.set_trace(context=10);

def simpler():
    #                    q0  q1  q2  q3
    HB = (tensor_product(IN,  I,  I,  I) +
          tensor_product( I, IN,  I,  I) +
          tensor_product( I,  I, IN,  I) +
          tensor_product( I,  I,  I, IN))

    #                    q0  q1  q2  q3
    HP = (tensor_product(P0, P0, P0,  I) +
          tensor_product( I, P1, P0, P0))

    evolve(HB, HP)

if __name__ == "__main__":
#     main()
#     other()
    simpler()