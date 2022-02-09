import sympy
from term import Term, ComposedTerm
from abstract_hamiltonian import AbstractHamiltonian


class DirectHamiltonian(AbstractHamiltonian):

    def __init__(self, N):
        super().__init__(N)

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
        return self.to_string(final)

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
        return self.to_string(expression)

if __name__ == "__main__":
    for N in range(11, 30, 4):
        H = DirectHamiltonian(N).with_parameters(letters=False, show_products=False)
        print(str(H))
        exp0 = H.from_expansion()
        print(f"H = {exp0}\n")
        if H.show_products:
            exp0 = str(sympy.simplify(exp0))
            print(f"H = {exp0}\n")

        # exp0 = H.from_function()
        # exp0 = sympy.expand(exp0)
        # exp0 = str(exp0).replace("**2", "")
        # exp0 = str(sympy.simplify(exp0))
        # exp1 = str(sympy.simplify(H.from_expansion()))
        # exp2 = str(sympy.simplify(H.from_formula()))
        # assert(exp0 == exp1 and exp0 == exp2)
        # print(exp2)
