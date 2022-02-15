from sympy import multiplicity
from term import ComposedTermWithReplacedVariables, Term, ComposedTerm
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

        x = self.to_string(xs)
        y = self.to_string(ys)
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

    def quadratized_by_ishikawa(self):
        N = self.N
        nx = self.nx
        ny = self.ny
        b = nx + ny
        expression = []

        for i in range(1, nx+1):
            for j in range(i+1, nx+1):
                for k in range(1, ny+1):
                    for l in range(k+1, ny+1):
                        b += 1
                        exp2 = (i+j+k+l+2)
                        terms = [Term({i, j},       exp2),
                                 Term({i, k+nx},    exp2),
                                 Term({i, l+nx},    exp2),
                                 Term({j, k+nx},    exp2),
                                 Term({j, l+nx},    exp2),
                                 Term({k+nx, l+nx}, exp2),
                                 Term({k+nx, l+nx}, exp2),
                                 Term({i, b},       exp2, -2),
                                 Term({j, b},       exp2, -2),
                                 Term({k+nx, b},    exp2, -2),
                                 Term({l+nx, b},    exp2, -2),
                                 Term({b},          exp2, +3)]
                        expression += terms
        for i in range(1, nx+1):
            for j in range(i+1, nx+1):
                for k in range(1, ny+1):
                    b += 1
                    temp1 = Term(two_exponent=(i+j+2*k+1))
                    temp2 = Term(two_exponent=(i+j+k+2))
                    terms1 = [ComposedTermWithReplacedVariables({i, j},    temp1.copy(), temp2.copy()),
                              ComposedTermWithReplacedVariables({i, k+nx}, temp1.copy(), temp2.copy()),
                              ComposedTermWithReplacedVariables({j, k+nx}, temp1.copy(), temp2.copy()),
                              ComposedTermWithReplacedVariables({b},       temp1.copy(), temp2.copy())]
                    temp1.multiply(-1)
                    temp2.multiply(-1)
                    terms2 = [ComposedTermWithReplacedVariables({i, b},    temp1.copy(), temp2.copy()),
                              ComposedTermWithReplacedVariables({j, b},    temp1.copy(), temp2.copy()),
                              ComposedTermWithReplacedVariables({k+nx, b}, temp1.copy(), temp2.copy())]
                    expression += terms1 + terms2
        for i in range(1, nx+1):
            for j in range(1, ny+1):
                for k in range(j+1, ny+1):
                    b += 1
                    temp1 = Term(two_exponent=(2*i+j+k+1))
                    temp2 = Term(two_exponent=(i+j+k+2))
                    terms1 = [ComposedTermWithReplacedVariables({i, j+nx},    temp1.copy(), temp2.copy()),
                              ComposedTermWithReplacedVariables({i, k+nx},    temp1.copy(), temp2.copy()),
                              ComposedTermWithReplacedVariables({j+nx, k+nx}, temp1.copy(), temp2.copy()),
                              ComposedTermWithReplacedVariables({b},          temp1.copy(), temp2.copy())]
                    temp1.multiply(-1)
                    temp2.multiply(-1)
                    terms2 = [ComposedTermWithReplacedVariables({i, b},       temp1.copy(), temp2.copy()),
                              ComposedTermWithReplacedVariables({j+nx, b},    temp1.copy(), temp2.copy()),
                              ComposedTermWithReplacedVariables({k+nx, b},    temp1.copy(), temp2.copy())]
                    expression += terms1 + terms2
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

    def quadratized_by_dattani_chau(self):
        N = self.N
        nx = self.nx
        ny = self.ny
        xxyy = []
        xxy = []
        xyy = []
        for i in range(1, nx+1):
            for j in range(i+1, nx+1):
                for k in range(1, ny+1):
                    for l in range(k+1, ny+1):
                        term = Term({i, j, k+nx, l+nx}, (i+j+k+l+2))
                        term.children = []
                        xxyy.append(term)
        for i in range(1, nx+1):
            for j in range(i+1, nx+1):
                for k in range(1, ny+1):
                    term1 = Term({i, j, k+nx}, (i+j+2*k+1))
                    term2 = Term({i, j, k+nx}, (i+j+k+2))
                    term = ComposedTerm(term1, term2)
                    xxy.append(term)
        for i in range(1, nx+1):
            for j in range(1, ny+1):
                for k in range(j+1, ny+1):
                    term1 = Term({i, j+nx, k+nx}, (2*i+j+k+1))
                    term2 = Term({i, j+nx, k+nx}, (i+j+k+2))
                    term = ComposedTerm(term1, term2)
                    xyy.append(term)

        for child in xxy:
            for parent in xxyy:
                if child.vars.issubset(parent.vars):
                    parent.children.append(child)
                    break
        for child in xyy:
            for parent in xxyy:
                if child.vars.issubset(parent.vars):
                    parent.children.append(child)
                    break

        expression = []
        for i in range(len(xxyy)):
            expression += xxyy[i].quadratize_by_first_lemma(var=(nx+ny+i+1))

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
