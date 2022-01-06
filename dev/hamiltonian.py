import sympy
import numpy as np
from utils.utils import number_of_bits_for_x, number_of_bits_for_y


class Term:

    def __init__(self, variable_set=set(), two_exponent=0, multiplier=1):
        self.vars = variable_set
        self.exp2 = two_exponent
        self.mult = multiplier

    def __mul__(self, other):
        variable_set = self.vars.union(other.vars)
        two_exponent = self.exp2 + other.exp2
        multiplier   = self.mult * other.mult
        return Term(variable_set, two_exponent, multiplier)

    def same_monomial(self, other):
        return self.vars == other.vars

    def same_exponent(self, other):
        return self.exp2 == other.exp2

    def __add__(self, other):
        if not (self.same_monomial(other) and self.same_exponent(other)):
            raise RuntimeError()
        variable_set = self.vars.copy()
        two_exponent = self.exp2
        multiplier   = self.mult + other.mult
        return Term(variable_set, two_exponent, multiplier)

    def __sub__(self, other):
        copy = other.copy()
        copy.mult *= -1
        return self.__add__(copy)

    def copy(self):
        return Term(self.vars.copy(), self.exp2, self.mult)

    def variables_string(self, nx=None, letters=False, latex=False, show_products=False):
        if letters:
            if show_products:
                return "*".join([chr(96+v) for v in self.vars])
            return "".join([chr(96+v) for v in self.vars])
        if nx is not None:
            if latex:
                return "".join([f"x_{{v}}" if v < nx+1 else f"y_{{v-nx}}" for v in self.vars])
            return "".join([f"x{v}" if v < nx+1 else f"y{v-nx}" for v in self.vars])
        if latex:
            return "".join([f"b_{{v}}" for v in self.vars])
        return "".join([f"b{v}" for v in self.vars])

    def coefficient_string(self, ignore_variables=False):
        if self.exp2 != 0:
            if self.mult != 1:
                return f"{self.mult:+d}*2^{self.exp2}"
            else:
                return f"+2^{self.exp2}"
        if ignore_variables:
            return f"{self.mult:+d}"
        else:
            if self.mult == 1:
                return "+"
        return f"{self.mult:+d}"

    def as_string(self, nx=None, letters=False, latex=False, show_products=False):
        if self.mult == 0:
            return "+0"
        vars = self.variables_string(nx, letters, latex, show_products)
        if vars:
            if show_products:
                return f"{self.coefficient_string()}*{vars}"
            return f"{self.coefficient_string()}{vars}"
        return f"{self.coefficient_string(True)}"

    def __str__(self):
        return self.as_string()

class ComposedTerm:

    def __init__(self, term, *terms):
        self.vars = term.vars
        self.terms = [term] + list(terms)
        self.assert_vars()

    def same_monomial(self, other):
        return self.vars == other.vars

    def assert_vars(self):
        for term in self.terms:
            assert(self.vars == term.vars)

    def append(self, term):
        assert(self.vars == term.vars)
        self.terms.append(term)

    def variables_string(self, nx=None, letters=False, latex=False, show_products=False):
        if letters:
            if show_products:
                return "*".join([chr(96+v) for v in self.vars])
            return "".join([chr(96+v) for v in self.vars])
        if nx is not None:
            if latex:
                return "".join([f"x_{{v}}" if v < nx+1 else f"y_{{v-nx}}" for v in self.vars])
            return "".join([f"x{v}" if v < nx+1 else f"y{v-nx}" for v in self.vars])
        if latex:
            return "".join([f"b_{{v}}" for v in self.vars])
        return "".join([f"b{v}" for v in self.vars])

    def as_string(self, nx=None, letters=False, latex=False, show_products=False):
        if len(self.terms) == 1:
            return self.terms[0].as_string(nx, letters, latex, show_products)
        coef = f"+({''.join(term.coefficient_string(True) for term in self.terms)})"
        vars = self.variables_string(nx, letters, latex, show_products)
        if vars:
            if show_products:
                return f"{coef}*{vars}"
            return f"{coef}{vars}"
        return f"{coef}"

    def __str__(self):
        return self.as_string()

class DirectHamiltonian:

    def __init__(self, N):
        self.N  = N
        self.n  = int(np.ceil(np.log2(N)))
        self.nx = int(number_of_bits_for_x(N))
        self.ny = int(number_of_bits_for_y(N))

    def from_function(self):
        N = self.N
        nx = self.nx
        ny = self.ny

        xs = [Term()]; ys = [Term()]
        for x in range(1, nx+1):
            xs.append(Term({x}, x))
        for y in range(1, ny+1):
            ys.append(Term({y+nx}, y))

        xs.reverse()
        ys.reverse()

        x = "".join([t.as_string(letters=True, show_products=True) for t in xs])
        y = "".join([t.as_string(letters=True, show_products=True) for t in ys])
        return f"({N} - ({x})*({y}))^2"

    def from_expansion(self):
        nx = self.nx
        ny = self.ny

        xs = []; ys = []
        for x in range(1, nx+1):
            xs.append(Term({x}, x))
        for y in range(1, ny+1):
            ys.append(Term({y+nx}, y))

        terms = []
        for x in xs:
            for y in ys:
                terms.append(x * y)
        terms += xs + ys
        for term in terms:
            term.mult = -1
        terms.append(Term(multiplier=(self.N-1)))

        expression = []
        for term1 in terms:
            for term2 in terms:
                expression.append(term1 * term2)

        expression.sort(key=lambda t: [len(t.vars)] + sorted(list(t.vars)))

        simplified = []
        term = expression[0]
        for i in range(1, len(expression)):
            if term.same_monomial(expression[i]) and term.same_exponent(expression[i]):
                term += expression[i]
            else:
                simplified.append(term)
                term = expression[i]
        simplified.append(term)

        final = []
        term = ComposedTerm(simplified[0])
        for i in range(1, len(simplified)):
            if term.same_monomial(simplified[i]):
                term.append(simplified[i])
            else:
                final.append(term)
                term = ComposedTerm(simplified[i])
        final.append(term)
        return self.expression_to_string(final)

    def from_formula(self):
        N = self.N
        nx = self.nx
        ny = self.ny
        expression = []
        for i in range(1, nx+1):
            for j in range(i+1, nx+1):
                for k in range(1, ny+1):
                    for l in range(k+1, ny+1):
                        term = Term({i, j, k+nx, l+nx}, (i+j+k+l+2))
                        expression.append(term)
        for i in range(1, nx+1):
            for j in range(i+1, nx+1):
                for k in range(1, ny+1):
                    term1 = Term({i, j, k+nx}, (i+j+2*k+1))
                    term2 = Term({i, j, k+nx}, (i+j+k+2))
                    term = ComposedTerm(term1, term2)
                    expression.append(term)
        for i in range(1, nx+1):
            for j in range(1, ny+1):
                for k in range(j+1, ny+1):
                    term1 = Term({i, j+nx, k+nx}, (2*i+j+k+1))
                    term2 = Term({i, j+nx, k+nx}, (i+j+k+2))
                    term = ComposedTerm(term1, term2)
                    expression.append(term)
        for i in range(1, nx+1):
            for j in range(1, ny+1):
                term1 = Term({i, j+nx}, (2*i+2*j))
                term2 = Term({i, j+nx}, (2*i+j+1))
                term3 = Term({i, j+nx}, (i+2*j+1))
                term4 = Term({i, j+nx}, (i+j+1), (2-N))
                term = ComposedTerm(term1, term2, term3, term4)
                expression.append(term)
        for i in range(1, nx+1):
            for j in range(i+1, nx+1):
                term = Term({i, j}, (i+j+1))
                expression.append(term)
        for i in range(1, ny+1):
            for j in range(i+1, ny+1):
                term = Term({i+nx, j+nx}, (i+j+1))
                expression.append(term)
        for i in range(1, nx+1):
            term1 = Term({i}, (2*i))
            term2 = Term({i}, (i+1), (1-N))
            term = ComposedTerm(term1, term2)
            expression.append(term)
        for i in range(1, ny+1):
            term1 = Term({i+nx}, (2*i))
            term2 = Term({i+nx}, (i+1), (1-N))
            term = ComposedTerm(term1, term2)
            expression.append(term)
        expression.append(Term(multiplier=(N-1)**2))
        expression.sort(key=lambda t: [len(t.vars)] + sorted(list(t.vars)))
        return self.expression_to_string(expression)

    def expression_to_string(self, expression):
        return "".join([t.as_string(letters=True, show_products=True) for t in expression])

    def as_string(self):
        return f"H: N={self.N}, n={self.n}, nx={self.nx}, ny={self.ny}"

    def __str__(self):
        return self.as_string()

if __name__ == "__main__":
    for N in range(11, 55, 2):
        hamiltonian = DirectHamiltonian(N)
        exp0 = hamiltonian.from_function()
        exp0 = sympy.expand(exp0)
        exp0 = str(exp0).replace("**2", "")
        exp0 = str(sympy.simplify(exp0))
        exp1 = str(sympy.simplify(hamiltonian.from_expansion()))
        exp2 = str(sympy.simplify(hamiltonian.from_formula()))
        assert(exp0 == exp1 and exp0 == exp2)
        print(str(hamiltonian))
        print(exp2)
