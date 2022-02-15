
class AbstractTerm:

    def __init__(self):
        self.children = []

    def same_monomial(self, other):
        return self.vars == other.vars

    def variables_string(self, nx=None, ny=None, letters=False, latex=False, show_products=False):
        join = "*" if show_products else ""
        if letters:
            return join.join([chr(96+v) for v in self.vars])
        if nx is not None and ny is not None:
            if latex:
                return join.join([f"x_{{v:02d}}" if v <= nx else (f"y_{{v-nx:02d}}" if v <= (nx+ny) else f"b_{{v-nx-ny:02d}}") for v in self.vars])
            return join.join([f"x{v:02d}" if v <= nx else (f"y{v-nx:02d}" if v <= (nx+ny) else f"b{v-nx-ny:02d}") for v in self.vars])
        if latex:
            return join.join([f"b_{{v:02d}}" for v in self.vars])
        return join.join([f"b{v:02d}" for v in self.vars])

    def __str__(self):
        return self.as_string()

class Term(AbstractTerm):

    def __init__(self, variable_set=set(), two_exponent=0, multiplier=1):
        super().__init__()
        self.vars = variable_set
        self.exp2 = two_exponent
        self.mult = multiplier

    def same_exponent(self, other):
        return self.exp2 == other.exp2

    def multiply(self, mult):
        self.mult *= mult

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

    def quadratize_by_first_lemma(self, var):
        copy = self.copy()
        copy.multiply(3)
        children = [child.copy() for child in self.children]
        term1 = [ComposedTermWithReplacedVariables({var}, copy, *children)]

        term2 = []
        term3 = []
        vars = list(self.vars)
        for i in range(len(vars)):
            for j in range(i+1, len(vars)):
                copy = self.copy()
                copy.vars = {vars[i], vars[j]}
                term2.append(copy)

                children = []
                for child in self.children:
                    if vars[i] in child.vars and vars[j] in child.vars:
                        children.append(child.copy())
                if children:
                    composed = ComposedTermWithReplacedVariables({vars[i], vars[j]}, *children)
                    term3.append(composed)

        term4 = []
        for i in range(len(vars)):
            copy = self.copy()
            copy.multiply(-2)

            children = []
            for child in self.children:
                if vars[i] in child.vars:
                    child_ = child.copy()
                    child_.multiply(-1)
                    children.append(child_)
            if children:
                composed = ComposedTermWithReplacedVariables({vars[i], var}, copy, *children)
                term4.append(composed)
        return term1 + term2 + term3 + term4

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

    def as_string(self, nx=None, ny=None, letters=False, latex=False, show_products=False):
        if self.mult == 0:
            return "+0"
        vars = self.variables_string(nx, ny, letters, latex, show_products)
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

    def multiply(self, mult):
        for term in self.terms:
            term.multiply(mult)

    def copy(self):
        terms = [term.copy() for term in self.terms]
        return ComposedTerm(*terms)

    def as_string(self, nx=None, ny=None, letters=False, latex=False, show_products=False):
        if len(self.terms) == 1:
            return self.terms[0].as_string(nx, ny, letters, latex, show_products)
        coef = f"+({''.join(term.coefficient_string(True) for term in self.terms)})"
        vars = self.variables_string(nx, ny, letters, latex, show_products)
        if vars:
            if show_products:
                return f"{coef}*{vars}"
            return f"{coef}{vars}"
        return f"{coef}"

class ComposedTermWithReplacedVariables(ComposedTerm):

    def __init__(self, vars, *terms):
        self.vars = vars
        self.terms = []
        for term in terms:
            if isinstance(term, ComposedTerm):
                self.terms += term.terms
            else:
                self.terms.append(term)

    def assert_vars(self):
        pass

    def append(self, term):
        self.terms.append(term)
