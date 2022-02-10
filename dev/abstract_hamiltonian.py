import numpy as np
from utils.utils import number_of_bits_for_x
from utils.utils import number_of_bits_for_y


class AbstractHamiltonian:

    def __init__(self, N):
        self.N  = N
        self.n  = int(np.ceil(np.log2(N)))
        self.nx = int(number_of_bits_for_x(N))
        self.ny = int(number_of_bits_for_y(N))

        self.letters = False
        self.latex = False
        self.show_products = False

    def as_string(self):
        return f"H: N={self.N}, n={self.n}, nx={self.nx}, ny={self.ny}"

    def __str__(self):
        return self.as_string()

    def to_string(self, expression):
        return "".join([t.as_string(self.nx, self.ny, self.letters, self.latex, self.show_products) for t in expression])

    def with_parameters(self, letters=False, latex=False, show_products=False):
        self.letters = letters
        self.latex = latex
        self.show_products = show_products
        return self
