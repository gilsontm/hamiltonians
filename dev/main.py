import sympy
from direct_hamiltonian_sympy import DirectHamiltonianSympy


def minimize(expression, H):
    exp = sympy.sympify(expression)
    symbols = list(exp.free_symbols)
    symbols.sort(key=lambda s: s.name, reverse=True)
    fun = sympy.lambdify(symbols, exp, "numpy")

    minimum = None
    for i in range(2**(len(symbols))):
        attr = bin(i)[2:].rjust(len(symbols), "0")
        energy = fun(*map(int, attr))
        if i == 0 or energy < minimum[0]:
            minimum = (energy, attr)

    y = int(minimum[1][:H.ny] + "1", base=2)
    x = int(minimum[1][H.ny:H.ny+H.nx] + "1", base=2)
    return symbols, minimum, x, y

def main():
    for N in range(81, 90, 2):
        H = DirectHamiltonianSympy(N).with_parameters(letters=False, show_products=True)
        print(H)
        exp0 = str(sympy.expand(H.from_simpy_formula()))
        exp1 = str(sympy.expand(H.from_alternate_formula()))
        sym1, res1, x1, y1 = minimize(exp1, H)
        assert(exp0 == exp1)

        # exp2 = str(sympy.expand(H.from_ishikawa_formula()))
        # print(exp2)
        # sym2, res2, x2, y2 = minimize(exp2, H)
        # print(f"[ISHI] x={x2:3d} y={y2:3d} E={res2[0]}")
        # assert(res1[0] == res2[0])
        # assert(x1 == x2)
        # assert(y1 == y2)

        exp3 = H.from_dattani_formula()
        exp3 = str(sympy.expand(exp3))
        print(exp3)
        sym3, res3, x3, y3 = minimize(exp3, H)
        print(f"[DATT] x={x3:3d} y={y3:3d} E={res3[0]}")
        print(sym3)

        assert(res1[0] == res3[0])
        assert(x1 == x3)
        assert(y1 == y3)

if __name__ == "__main__":
    main()
