from abstract_hamiltonian import AbstractHamiltonian


class DirectHamiltonianSympy(AbstractHamiltonian):

    def __init__(self, N):
        super().__init__(N)

    def from_simpy_formula(self):
        N = self.N
        nx = self.nx
        ny = self.ny
        expression = []
        for i in range(1, nx+1):
            for j in range(i+1, nx+1):
                for k in range(1, ny+1):
                    for l in range(k+1, ny+1):
                        term = f"+2^{i+j+k+l+2}*x{i:02d}*x{j:02d}*y{k:02d}*y{l:02d}"
                        expression.append(term)
        for i in range(1, nx+1):
            for j in range(i+1, nx+1):
                for k in range(1, ny+1):
                    term = f"+(2^{i+j+2*k+1}+2^{i+j+k+2})*x{i:02d}*x{j:02d}*y{k:02d}"
                    expression.append(term)
        for i in range(1, nx+1):
            for k in range(1, ny+1):
                for l in range(k+1, ny+1):
                    term = f"+(2^{2*i+k+l+1}+2^{i+k+l+2})*x{i:02d}*y{k:02d}*y{l:02d}"
                    expression.append(term)
        for i in range(1, nx+1):
            for k in range(1, ny+1):
                term = f"+(2^{2*i+2*k}+2^{2*i+k+1}+2^{i+2*k+1}-2^{i+k+1}*{N-2})*x{i:02d}*y{k:02d}"
                expression.append(term)
        for i in range(1, nx+1):
            for j in range(i+1, nx+1):
                term = f"+2^{i+j+1}*x{i:02d}*x{j:02d}"
                expression.append(term)
        for k in range(1, ny+1):
            for l in range(k+1, ny+1):
                term = f"+2^{k+l+1}*y{k:02d}*y{l:02d}"
                expression.append(term)
        for i in range(1, nx+1):
            term = f"+(2^{2*i}-2^{i+1}*{N-1})*x{i:02d}"
            expression.append(term)
        for k in range(1, ny+1):
            term = f"+(2^{2*k}-2^{k+1}*{N-1})*y{k:02d}"
            expression.append(term)
        expression.append(f"+{N-1}^2")
        return "".join(expression)

    def from_alternate_formula(self):
        N = self.N
        nx = self.nx
        ny = self.ny

        term1 = f"+2^{1+2+1+2+2}*x{1:02d}*x{2:02d}*y{1:02d}*y{2:02d}"
        term2 = f"+(2^{1+2+2*1+1}+2^{1+2+1+2})*x{1:02d}*x{2:02d}*y{1:02d}"
        term3 = f"+(2^{1+2+2*2+1}+2^{1+2+2+2})*x{1:02d}*x{2:02d}*y{2:02d}"
        term4 = f"+(2^{2*1+1+2+1}+2^{1+1+2+2})*x{1:02d}*y{1:02d}*y{2:02d}"
        term5 = f"+(2^{2*2+1+2+1}+2^{2+1+2+2})*x{2:02d}*y{1:02d}*y{2:02d}"
        expression = [term1, term2, term3, term4, term5]

        for l in range(3, ny+1):
            term1 = f"+2^{1+2+1+l+2}*x{1:02d}*x{2:02d}*y{1:02d}*y{l:02d}"
            term2 = f"+(2^{1+2+2*l+1}+2^{1+2+l+2})*x{1:02d}*x{2:02d}*y{l:02d}"
            term3 = f"+(2^{2*1+1+l+1}+2^{1+1+l+2})*x{1:02d}*y{1:02d}*y{l:02d}"
            term4 = f"+(2^{2*2+1+l+1}+2^{2+1+l+2})*x{2:02d}*y{1:02d}*y{l:02d}"
            expression += [term1, term2, term3, term4]
        for j in range(3, nx+1):
            term1 = f"+2^{1+j+1+2+2}*x{1:02d}*x{j:02d}*y{1:02d}*y{2:02d}"
            term2 = f"+(2^{2*j+1+2+1}+2^{j+1+2+2})*x{j:02d}*y{1:02d}*y{2:02d}"
            term3 = f"+(2^{1+j+2*1+1}+2^{1+j+1+2})*x{1:02d}*x{j:02d}*y{1:02d}"
            term4 = f"+(2^{1+j+2*2+1}+2^{1+j+2+2})*x{1:02d}*x{j:02d}*y{2:02d}"
            expression += [term1, term2, term3, term4]

        for j in range(3, nx+1):
            for l in range(3, ny+1):
                term1 = f"+2^{1+j+1+l+2}*x{1:02d}*x{j:02d}*y{1:02d}*y{l:02d}"
                term2 = f"+(2^{2*j+1+l+1}+2^{j+1+l+2})*x{j:02d}*y{1:02d}*y{l:02d}"
                term3 = f"+(2^{1+j+2*l+1}+2^{1+j+l+2})*x{1:02d}*x{j:02d}*y{l:02d}"
                expression += [term1, term2, term3]

        for k in range(2, ny+1):
            for l in range(k+1, ny+1):
                term1 = f"+2^{1+2+k+l+2}*x{1:02d}*x{2:02d}*y{k:02d}*y{l:02d}"
                term2 = f"+(2^{2*1+k+l+1}+2^{1+k+l+2})*x{1:02d}*y{k:02d}*y{l:02d}"
                term3 = f"+(2^{2*2+k+l+1}+2^{2+k+l+2})*x{2:02d}*y{k:02d}*y{l:02d}"
                expression += [term1, term2, term3]
        for j in range(3, nx+1):
            for k in range(2, ny+1):
                for l in range(k+1, ny+1):
                    term1 = f"+2^{1+j+k+l+2}*x{1:02d}*x{j:02d}*y{k:02d}*y{l:02d}"
                    term2 = f"+(2^{2*j+k+l+1}+2^{j+k+l+2})*x{j:02d}*y{k:02d}*y{l:02d}"
                    expression += [term1, term2]
        for i in range(2, nx+1):
            for j in range(i+1, nx+1):
                term1 = f"+2^{i+j+1+2+2}*x{i:02d}*x{j:02d}*y{1:02d}*y{2:02d}"
                term2 = f"+(2^{i+j+2*1+1}+2^{i+j+1+2})*x{i:02d}*x{j:02d}*y{1:02d}"
                term3 = f"+(2^{i+j+2*2+1}+2^{i+j+2+2})*x{i:02d}*x{j:02d}*y{2:02d}"
                expression += [term1, term2, term3]
        for i in range(2, nx+1):
            for j in range(i+1, nx+1):
                for l in range(3, ny+1):
                    term1 = f"+2^{i+j+1+l+2}*x{i:02d}*x{j:02d}*y{1:02d}*y{l:02d}"
                    term2 = f"+(2^{i+j+2*l+1}+2^{i+j+l+2})*x{i:02d}*x{j:02d}*y{l:02d}"
                    expression += [term1, term2]
        for i in range(2, nx+1):
            for j in range(i+1, nx+1):
                for k in range(2, ny+1):
                    for l in range(k+1, ny+1):
                        term = f"+2^{i+j+k+l+2}*x{i:02d}*x{j:02d}*y{k:02d}*y{l:02d}"
                        expression.append(term)

        for i in range(1, nx+1):
            for k in range(1, ny+1):
                term = f"+(2^{2*i+2*k}+2^{2*i+k+1}+2^{i+2*k+1}-2^{i+k+1}*{N-2})*x{i:02d}*y{k:02d}"
                expression.append(term)
        for i in range(1, nx+1):
            for j in range(i+1, nx+1):
                term = f"+2^{i+j+1}*x{i:02d}*x{j:02d}"
                expression.append(term)
        for k in range(1, ny+1):
            for l in range(k+1, ny+1):
                term = f"+2^{k+l+1}*y{k:02d}*y{l:02d}"
                expression.append(term)
        for i in range(1, nx+1):
            term = f"+(2^{2*i}-2^{i+1}*{N-1})*x{i:02d}"
            expression.append(term)
        for k in range(1, ny+1):
            term = f"+(2^{2*k}-2^{k+1}*{N-1})*y{k:02d}"
            expression.append(term)
        expression.append(f"+{N-1}^2")
        return "".join(expression)

    def from_ishikawa_formula(self):
        N = self.N
        nx = self.nx
        ny = self.ny
        expression = []
        m = 0
        for i in range(1, nx+1):
            for j in range(i+1, nx+1):
                for k in range(1, ny+1):
                    for l in range(k+1, ny+1):
                        # term = f"+2^{i+j+k+l+2}*x{i:02d}*x{j:02d}*y{k:02d}*y{l:02d}"
                        m += 1
                        term = f"+2^{i+j+k+l+2}*(x{i:02d}*x{j:02d} + x{i:02d}*y{k:02d} + x{i:02d}*y{l:02d} + x{j:02d}*y{k:02d} + x{j:02d}*y{l:02d} + y{k:02d}*y{l:02d} + b{m:02d}*(3-2*x{i:02d}-2*x{j:02d}-2*y{k:02d}-2*y{l:02d}))"
                        expression.append(term)
        for i in range(1, nx+1):
            for j in range(i+1, nx+1):
                for k in range(1, ny+1):
                    # term = f"+(2^{i+j+2*k+1}+2^{i+j+k+2})*x{i:02d}*x{j:02d}*y{k:02d}"
                    m += 1
                    term = f"+(2^{i+j+2*k+1}+2^{i+j+k+2})*(x{i:02d}*x{j:02d} + x{i:02d}*y{k:02d} + x{j:02d}*y{k:02d} + b{m:02d}*(1-x{i:02d}-x{j:02d}-y{k:02d}))"
                    expression.append(term)
        for i in range(1, nx+1):
            for k in range(1, ny+1):
                for l in range(k+1, ny+1):
                    # term = f"+(2^{2*i+k+l+1}+2^{i+k+l+2})*x{i:02d}*y{k:02d}*y{l:02d}"
                    m += 1
                    term = f"+(2^{2*i+k+l+1}+2^{i+k+l+2})*(x{i:02d}*y{k:02d} + x{i:02d}*y{l:02d} + y{k:02d}*y{l:02d} + b{m:02d}*(1-x{i:02d}-y{k:02d}-y{l:02d}))"
                    expression.append(term)
        for i in range(1, nx+1):
            for k in range(1, ny+1):
                term = f"+(2^{2*i+2*k}+2^{2*i+k+1}+2^{i+2*k+1}-2^{i+k+1}*{N-2})*x{i:02d}*y{k:02d}"
                expression.append(term)
        for i in range(1, nx+1):
            for j in range(i+1, nx+1):
                term = f"+2^{i+j+1}*x{i:02d}*x{j:02d}"
                expression.append(term)
        for k in range(1, ny+1):
            for l in range(k+1, ny+1):
                term = f"+2^{k+l+1}*y{k:02d}*y{l:02d}"
                expression.append(term)
        for i in range(1, nx+1):
            term = f"+(2^{2*i}-2^{i+1}*{N-1})*x{i:02d}"
            expression.append(term)
        for k in range(1, ny+1):
            term = f"+(2^{2*k}-2^{k+1}*{N-1})*y{k:02d}"
            expression.append(term)
        expression.append(f"+{N-1}^2")
        return "".join(expression)

    def from_dattani_formula(self):
        N = self.N
        nx = self.nx
        ny = self.ny

        assert(nx >= 2)
        assert(ny >= 2)

        m = 0
        m += 1
        term1 = f"+(3*2^{1+2+1+2+2}+(2^{1+2+2*1+1}+2^{1+2+1+2})+(2^{1+2+2*2+1}+2^{1+2+2+2})+(2^{2*1+1+2+1}+2^{1+1+2+2})+(2^{2*2+1+2+1}+2^{2+1+2+2}))*b{m:02d}"

        term2 = f"+(2^{1+2+1+2+2}+(2^{1+2+2*1+1}+2^{1+2+1+2})+(2^{1+2+2*2+1}+2^{1+2+2+2}))*x{1:02d}*x{2:02d}"
        term3 = f"+(2^{1+2+1+2+2}+(2^{1+2+2*1+1}+2^{1+2+1+2})+(2^{2*1+1+2+1}+2^{1+1+2+2}))*x{1:02d}*y{1:02d}"
        term4 = f"+(2^{1+2+1+2+2}+(2^{1+2+2*2+1}+2^{1+2+2+2})+(2^{2*1+1+2+1}+2^{1+1+2+2}))*x{1:02d}*y{2:02d}"
        term5 = f"+(2^{1+2+1+2+2}+(2^{1+2+2*1+1}+2^{1+2+1+2})+(2^{2*2+1+2+1}+2^{2+1+2+2}))*x{2:02d}*y{1:02d}"
        term6 = f"+(2^{1+2+1+2+2}+(2^{1+2+2*2+1}+2^{1+2+2+2})+(2^{2*2+1+2+1}+2^{2+1+2+2}))*x{2:02d}*y{2:02d}"
        term7 = f"+(2^{1+2+1+2+2}+(2^{2*1+1+2+1}+2^{1+1+2+2})+(2^{2*2+1+2+1}+2^{2+1+2+2}))*y{1:02d}*y{2:02d}"

        term8  = f"-(2*2^{1+2+1+2+2}+(2^{1+2+2*1+1}+2^{1+2+1+2})+(2^{1+2+2*2+1}+2^{1+2+2+2})+(2^{2*1+1+2+1}+2^{1+1+2+2}))*x{1:02d}*b{m:02d}"
        term9  = f"-(2*2^{1+2+1+2+2}+(2^{1+2+2*1+1}+2^{1+2+1+2})+(2^{1+2+2*2+1}+2^{1+2+2+2})+(2^{2*2+1+2+1}+2^{2+1+2+2}))*x{2:02d}*b{m:02d}"
        term10 = f"-(2*2^{1+2+1+2+2}+(2^{1+2+2*1+1}+2^{1+2+1+2})+(2^{2*1+1+2+1}+2^{1+1+2+2})+(2^{2*2+1+2+1}+2^{2+1+2+2}))*y{1:02d}*b{m:02d}"
        term11 = f"-(2*2^{1+2+1+2+2}+(2^{1+2+2*2+1}+2^{1+2+2+2})+(2^{2*1+1+2+1}+2^{1+1+2+2})+(2^{2*2+1+2+1}+2^{2+1+2+2}))*y{2:02d}*b{m:02d}"
        expression = [term1, term2, term3, term4, term5, term6, term7, term8, term9, term10, term11]

        for l in range(3, ny+1):
            m += 1
            term1 = f"+(3*2^{1+2+1+l+2}+(2^{1+2+2*l+1}+2^{1+2+l+2})+(2^{2*1+1+l+1}+2^{1+1+l+2})+(2^{2*2+1+l+1}+2^{2+1+l+2}))*b{m:02d}"

            term2 = f"+(2^{1+2+1+l+2}+(2^{1+2+2*l+1}+2^{1+2+l+2}))*x{1:02d}*x{2:02d}"
            term3 = f"+(2^{1+2+1+l+2}+(2^{2*1+1+l+1}+2^{1+1+l+2}))*x{1:02d}*y{1:02d}"
            term4 = f"+(2^{1+2+1+l+2}+(2^{1+2+2*l+1}+2^{1+2+l+2})+(2^{2*1+1+l+1}+2^{1+1+l+2}))*x{1:02d}*y{l:02d}"
            term5 = f"+(2^{1+2+1+l+2}+(2^{2*2+1+l+1}+2^{2+1+l+2}))*x{2:02d}*y{1:02d}"
            term6 = f"+(2^{1+2+1+l+2}+(2^{1+2+2*l+1}+2^{1+2+l+2})+(2^{2*2+1+l+1}+2^{2+1+l+2}))*x{2:02d}*y{l:02d}"
            term7 = f"+(2^{1+2+1+l+2}+(2^{2*1+1+l+1}+2^{1+1+l+2})+(2^{2*2+1+l+1}+2^{2+1+l+2}))*y{1:02d}*y{l:02d}"

            term8  = f"-(2*2^{1+2+1+l+2}+(2^{1+2+2*l+1}+2^{1+2+l+2})+(2^{2*1+1+l+1}+2^{1+1+l+2}))*x{1:02d}*b{m:02d}"
            term9  = f"-(2*2^{1+2+1+l+2}+(2^{1+2+2*l+1}+2^{1+2+l+2})+(2^{2*2+1+l+1}+2^{2+1+l+2}))*x{2:02d}*b{m:02d}"
            term10 = f"-(2*2^{1+2+1+l+2}+(2^{2*1+1+l+1}+2^{1+1+l+2})+(2^{2*2+1+l+1}+2^{2+1+l+2}))*y{1:02d}*b{m:02d}"
            term11 = f"-(2*2^{1+2+1+l+2}+(2^{1+2+2*l+1}+2^{1+2+l+2})+(2^{2*1+1+l+1}+2^{1+1+l+2})+(2^{2*2+1+l+1}+2^{2+1+l+2}))*y{l:02d}*b{m:02d}"
            expression += [term1, term2, term3, term4, term5, term6, term7, term8, term9, term10, term11]
        for j in range(3, nx+1):
            m += 1
            term1 = f"+(3*2^{1+j+1+2+2}+(2^{2*j+1+2+1}+2^{j+1+2+2})+(2^{1+j+2*1+1}+2^{1+j+1+2})+(2^{1+j+2*2+1}+2^{1+j+2+2}))*b{m:02d}"

            term2 = f"+(2^{1+j+1+2+2}+(2^{1+j+2*1+1}+2^{1+j+1+2})+(2^{1+j+2*2+1}+2^{1+j+2+2}))*x{1:02d}*x{j:02d}"
            term3 = f"+(2^{1+j+1+2+2}+(2^{1+j+2*1+1}+2^{1+j+1+2}))*x{1:02d}*y{1:02d}"
            term4 = f"+(2^{1+j+1+2+2}+(2^{1+j+2*2+1}+2^{1+j+2+2}))*x{1:02d}*y{2:02d}"
            term5 = f"+(2^{1+j+1+2+2}+(2^{2*j+1+2+1}+2^{j+1+2+2})+(2^{1+j+2*1+1}+2^{1+j+1+2}))*x{j:02d}*y{1:02d}"
            term6 = f"+(2^{1+j+1+2+2}+(2^{2*j+1+2+1}+2^{j+1+2+2})+(2^{1+j+2*2+1}+2^{1+j+2+2}))*x{j:02d}*y{2:02d}"
            term7 = f"+(2^{1+j+1+2+2}+(2^{2*j+1+2+1}+2^{j+1+2+2}))*y{1:02d}*y{2:02d}"

            term8  = f"-(2*2^{1+j+1+2+2}+(2^{1+j+2*1+1}+2^{1+j+1+2})+(2^{1+j+2*2+1}+2^{1+j+2+2}))*x{1:02d}*b{m:02d}"
            term9  = f"-(2*2^{1+j+1+2+2}+(2^{2*j+1+2+1}+2^{j+1+2+2})+(2^{1+j+2*1+1}+2^{1+j+1+2})+(2^{1+j+2*2+1}+2^{1+j+2+2}))*x{j:02d}*b{m:02d}"
            term10 = f"-(2*2^{1+j+1+2+2}+(2^{2*j+1+2+1}+2^{j+1+2+2})+(2^{1+j+2*1+1}+2^{1+j+1+2}))*y{1:02d}*b{m:02d}"
            term11 = f"-(2*2^{1+j+1+2+2}+(2^{2*j+1+2+1}+2^{j+1+2+2})+(2^{1+j+2*2+1}+2^{1+j+2+2}))*y{2:02d}*b{m:02d}"
            expression += [term1, term2, term3, term4, term5, term6, term7, term8, term9, term10, term11]

        for j in range(3, nx+1):
            for l in range(3, ny+1):
                m += 1
                term1 = f"+(3*2^{1+j+1+l+2}+(2^{2*j+1+l+1}+2^{j+1+l+2})+(2^{1+j+2*l+1}+2^{1+j+l+2}))*b{m:02d}"

                term2 = f"+(2^{1+j+1+l+2}+(2^{1+j+2*l+1}+2^{1+j+l+2}))*x{1:02d}*x{j:02d}"
                term3 = f"+(2^{1+j+1+l+2})*x{1:02d}*y{1:02d}"
                term4 = f"+(2^{1+j+1+l+2}+(2^{1+j+2*l+1}+2^{1+j+l+2}))*x{1:02d}*y{l:02d}"
                term5 = f"+(2^{1+j+1+l+2}+(2^{2*j+1+l+1}+2^{j+1+l+2}))*x{j:02d}*y{1:02d}"
                term6 = f"+(2^{1+j+1+l+2}+(2^{2*j+1+l+1}+2^{j+1+l+2})+(2^{1+j+2*l+1}+2^{1+j+l+2}))*x{j:02d}*y{l:02d}"
                term7 = f"+(2^{1+j+1+l+2}+(2^{2*j+1+l+1}+2^{j+1+l+2}))*y{1:02d}*y{l:02d}"

                term8  = f"-(2*2^{1+j+1+l+2}+(2^{1+j+2*l+1}+2^{1+j+l+2}))*x{1:02d}*b{m:02d}"
                term9  = f"-(2*2^{1+j+1+l+2}+(2^{2*j+1+l+1}+2^{j+1+l+2})+(2^{1+j+2*l+1}+2^{1+j+l+2}))*x{j:02d}*b{m:02d}"
                term10 = f"-(2*2^{1+j+1+l+2}+(2^{2*j+1+l+1}+2^{j+1+l+2}))*y{1:02d}*b{m:02d}"
                term11 = f"-(2*2^{1+j+1+l+2}+(2^{2*j+1+l+1}+2^{j+1+l+2})+(2^{1+j+2*l+1}+2^{1+j+l+2}))*y{l:02d}*b{m:02d}"
                expression += [term1, term2, term3, term4, term5, term6, term7, term8, term9, term10, term11]

        for k in range(2, ny+1):
            for l in range(k+1, ny+1):
                m += 1
                term1 = f"+(3*2^{1+2+k+l+2}+(2^{2*1+k+l+1}+2^{1+k+l+2})+(2^{2*2+k+l+1}+2^{2+k+l+2}))*b{m:02d}"

                term2 = f"+(2^{1+2+k+l+2})*x{1:02d}*x{2:02d}"
                term3 = f"+(2^{1+2+k+l+2}+(2^{2*1+k+l+1}+2^{1+k+l+2}))*x{1:02d}*y{k:02d}"
                term4 = f"+(2^{1+2+k+l+2}+(2^{2*1+k+l+1}+2^{1+k+l+2}))*x{1:02d}*y{l:02d}"
                term5 = f"+(2^{1+2+k+l+2}+(2^{2*2+k+l+1}+2^{2+k+l+2}))*x{2:02d}*y{k:02d}"
                term6 = f"+(2^{1+2+k+l+2}+(2^{2*2+k+l+1}+2^{2+k+l+2}))*x{2:02d}*y{l:02d}"
                term7 = f"+(2^{1+2+k+l+2}+(2^{2*1+k+l+1}+2^{1+k+l+2})+(2^{2*2+k+l+1}+2^{2+k+l+2}))*y{k:02d}*y{l:02d}"

                term8  = f"-(2*2^{1+2+k+l+2}+(2^{2*1+k+l+1}+2^{1+k+l+2}))*x{1:02d}*b{m:02d}"
                term9  = f"-(2*2^{1+2+k+l+2}+(2^{2*2+k+l+1}+2^{2+k+l+2}))*x{2:02d}*b{m:02d}"
                term10 = f"-(2*2^{1+2+k+l+2}+(2^{2*1+k+l+1}+2^{1+k+l+2})+(2^{2*2+k+l+1}+2^{2+k+l+2}))*y{k:02d}*b{m:02d}"
                term11 = f"-(2*2^{1+2+k+l+2}+(2^{2*1+k+l+1}+2^{1+k+l+2})+(2^{2*2+k+l+1}+2^{2+k+l+2}))*y{l:02d}*b{m:02d}"
                expression += [term1, term2, term3, term4, term5, term6, term7, term8, term9, term10, term11]
        for j in range(3, nx+1):
            for k in range(2, ny+1):
                for l in range(k+1, ny+1):
                    m += 1
                    term1 = f"+(3*2^{1+j+k+l+2}+(2^{2*j+k+l+1}+2^{j+k+l+2}))*b{m:02d}"

                    term2 = f"+(2^{1+j+k+l+2})*x{1:02d}*x{j:02d}"
                    term3 = f"+(2^{1+j+k+l+2})*x{1:02d}*y{k:02d}"
                    term4 = f"+(2^{1+j+k+l+2})*x{1:02d}*y{l:02d}"
                    term5 = f"+(2^{1+j+k+l+2}+(2^{2*j+k+l+1}+2^{j+k+l+2}))*x{j:02d}*y{k:02d}"
                    term6 = f"+(2^{1+j+k+l+2}+(2^{2*j+k+l+1}+2^{j+k+l+2}))*x{j:02d}*y{l:02d}"
                    term7 = f"+(2^{1+j+k+l+2}+(2^{2*j+k+l+1}+2^{j+k+l+2}))*y{k:02d}*y{l:02d}"

                    term8  = f"-(2*2^{1+j+k+l+2})*x{1:02d}*b{m:02d}"
                    term9  = f"-(2*2^{1+j+k+l+2}+(2^{2*j+k+l+1}+2^{j+k+l+2}))*x{j:02d}*b{m:02d}"
                    term10 = f"-(2*2^{1+j+k+l+2}+(2^{2*j+k+l+1}+2^{j+k+l+2}))*y{k:02d}*b{m:02d}"
                    term11 = f"-(2*2^{1+j+k+l+2}+(2^{2*j+k+l+1}+2^{j+k+l+2}))*y{l:02d}*b{m:02d}"
                    expression += [term1, term2, term3, term4, term5, term6, term7, term8, term9, term10, term11]

        for i in range(2, nx+1):
            for j in range(i+1, nx+1):
                m += 1
                term1 = f"+(3*2^{i+j+1+2+2}+(2^{i+j+2*1+1}+2^{i+j+1+2})+(2^{i+j+2*2+1}+2^{i+j+2+2}))*b{m:02d}"

                term2 = f"+(2^{i+j+1+2+2}+(2^{i+j+2*1+1}+2^{i+j+1+2})+(2^{i+j+2*2+1}+2^{i+j+2+2}))*x{i:02d}*x{j:02d}"
                term3 = f"+(2^{i+j+1+2+2}+(2^{i+j+2*1+1}+2^{i+j+1+2}))*x{i:02d}*y{1:02d}"
                term4 = f"+(2^{i+j+1+2+2}+(2^{i+j+2*2+1}+2^{i+j+2+2}))*x{i:02d}*y{2:02d}"
                term5 = f"+(2^{i+j+1+2+2}+(2^{i+j+2*1+1}+2^{i+j+1+2}))*x{j:02d}*y{1:02d}"
                term6 = f"+(2^{i+j+1+2+2}+(2^{i+j+2*2+1}+2^{i+j+2+2}))*x{j:02d}*y{2:02d}"
                term7 = f"+(2^{i+j+1+2+2})*y{1:02d}*y{2:02d}"

                term8  = f"-(2*2^{i+j+1+2+2}+(2^{i+j+2*1+1}+2^{i+j+1+2})+(2^{i+j+2*2+1}+2^{i+j+2+2}))*x{i:02d}*b{m:02d}"
                term9  = f"-(2*2^{i+j+1+2+2}+(2^{i+j+2*1+1}+2^{i+j+1+2})+(2^{i+j+2*2+1}+2^{i+j+2+2}))*x{j:02d}*b{m:02d}"
                term10 = f"-(2*2^{i+j+1+2+2}+(2^{i+j+2*1+1}+2^{i+j+1+2}))*y{1:02d}*b{m:02d}"
                term11 = f"-(2*2^{i+j+1+2+2}+(2^{i+j+2*2+1}+2^{i+j+2+2}))*y{2:02d}*b{m:02d}"
                expression += [term1, term2, term3, term4, term5, term6, term7, term8, term9, term10, term11]
        for i in range(2, nx+1):
            for j in range(i+1, nx+1):
                for l in range(3, ny+1):
                    m += 1
                    term1 = f"+(3*2^{i+j+1+l+2}+(2^{i+j+2*l+1}+2^{i+j+l+2}))*b{m:02d}"

                    term2 = f"+(2^{i+j+1+l+2}+(2^{i+j+2*l+1}+2^{i+j+l+2}))*x{i:02d}*x{j:02d}"
                    term3 = f"+(2^{i+j+1+l+2})*x{i:02d}*y{1:02d}"
                    term4 = f"+(2^{i+j+1+l+2}+(2^{i+j+2*l+1}+2^{i+j+l+2}))*x{i:02d}*y{l:02d}"
                    term5 = f"+(2^{i+j+1+l+2})*x{j:02d}*y{1:02d}"
                    term6 = f"+(2^{i+j+1+l+2}+(2^{i+j+2*l+1}+2^{i+j+l+2}))*x{j:02d}*y{l:02d}"
                    term7 = f"+(2^{i+j+1+l+2})*y{1:02d}*y{l:02d}"

                    term8  = f"-(2*2^{i+j+1+l+2}+(2^{i+j+2*l+1}+2^{i+j+l+2}))*x{i:02d}*b{m:02d}"
                    term9  = f"-(2*2^{i+j+1+l+2}+(2^{i+j+2*l+1}+2^{i+j+l+2}))*x{j:02d}*b{m:02d}"
                    term10 = f"-(2*2^{i+j+1+l+2})*y{1:02d}*b{m:02d}"
                    term11 = f"-(2*2^{i+j+1+l+2}+(2^{i+j+2*l+1}+2^{i+j+l+2}))*y{l:02d}*b{m:02d}"
                    expression += [term1, term2, term3, term4, term5, term6, term7, term8, term9, term10, term11]

        for i in range(2, nx+1):
            for j in range(i+1, nx+1):
                for k in range(2, ny+1):
                    for l in range(k+1, ny+1):
                        m += 1
                        term1 = f"+(3*2^{i+j+k+l+2})*b{m:02d}"

                        term2 = f"+(2^{i+j+k+l+2})*x{i:02d}*x{j:02d}"
                        term3 = f"+(2^{i+j+k+l+2})*x{i:02d}*y{k:02d}"
                        term4 = f"+(2^{i+j+k+l+2})*x{i:02d}*y{l:02d}"
                        term5 = f"+(2^{i+j+k+l+2})*x{j:02d}*y{k:02d}"
                        term6 = f"+(2^{i+j+k+l+2})*x{j:02d}*y{l:02d}"
                        term7 = f"+(2^{i+j+k+l+2})*y{k:02d}*y{l:02d}"

                        term8  = f"-(2*2^{i+j+k+l+2})*x{i:02d}*b{m:02d}"
                        term9  = f"-(2*2^{i+j+k+l+2})*x{j:02d}*b{m:02d}"
                        term10 = f"-(2*2^{i+j+k+l+2})*y{k:02d}*b{m:02d}"
                        term11 = f"-(2*2^{i+j+k+l+2})*y{l:02d}*b{m:02d}"
                        expression += [term1, term2, term3, term4, term5, term6, term7, term8, term9, term10, term11]

        for i in range(1, nx+1):
            for k in range(1, ny+1):
                term = f"+(2^{2*i+2*k}+2^{2*i+k+1}+2^{i+2*k+1}-2^{i+k+1}*{N-2})*x{i:02d}*y{k:02d}"
                expression.append(term)
        for i in range(1, nx+1):
            for j in range(i+1, nx+1):
                term = f"+2^{i+j+1}*x{i:02d}*x{j:02d}"
                expression.append(term)
        for k in range(1, ny+1):
            for l in range(k+1, ny+1):
                term = f"+2^{k+l+1}*y{k:02d}*y{l:02d}"
                expression.append(term)
        for i in range(1, nx+1):
            term = f"+(2^{2*i}-2^{i+1}*{N-1})*x{i:02d}"
            expression.append(term)
        for k in range(1, ny+1):
            term = f"+(2^{2*k}-2^{k+1}*{N-1})*y{k:02d}"
            expression.append(term)
        expression.append(f"+{N-1}^2")
        return "".join(expression)
