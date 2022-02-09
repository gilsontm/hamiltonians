
class AbstractTerm:

    def same_monomial(self, other):
        return self.vars == other.vars

    def variables_string(self, nx=None, letters=False, latex=False, show_products=False):
        join = "*" if show_products else ""
        if letters:
            return join.join([chr(96+v) for v in self.vars])
        if nx is not None:
            if latex:
                return join.join([f"x_{{v}}" if v < nx+1 else f"y_{{v-nx}}" for v in self.vars])
            return join.join([f"x{v}" if v < nx+1 else f"y{v-nx}" for v in self.vars])
        if latex:
            return join.join([f"b_{{v}}" for v in self.vars])
        return join.join([f"b{v}" for v in self.vars])

class Term(AbstractTerm):

    def __init__(self, variable_set=set(), two_exponent=0, multiplier=1):
        self.vars = variable_set
        self.exp2 = two_exponent
        self.mult = multiplier

    def same_exponent(self, other):
        return self.exp2 == other.exp2

    def __mul__(self, other):
        variable_set = self.vars.union(other.vars)
        two_exponent = self.exp2 + other.exp2
        multiplier   = self.mult * other.mult
        return Term(variable_set, two_exponent, multiplier)

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
                return "+1"
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

class ComposedTerm(AbstractTerm):

    def __init__(self, term, *terms):
        self.vars = term.vars
        self.terms = [term] + list(terms)
        self.assert_vars()

    def assert_vars(self):
        for term in self.terms:
            assert(self.vars == term.vars)

    def append(self, term):
        assert(self.vars == term.vars)
        self.terms.append(term)

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
