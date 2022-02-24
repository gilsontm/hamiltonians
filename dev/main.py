import sympy
import itertools
import numpy as np
from direct_hamiltonian import DirectHamiltonian
from direct_hamiltonian_sympy import DirectHamiltonianSympy


def minimize(expression, H):
    exp = sympy.sympify(expression)
    symbols = list(exp.free_symbols)
    symbols.sort(key=lambda s: s.name, reverse=True)
    fun = sympy.lambdify((symbols,), exp, "numpy")

    energy = 1e10
    window = 256
    iterator = itertools.product([0, 1], repeat=len(symbols))

    while True:
        combinations = np.bool_(list(itertools.islice(iterator, window))).transpose()
        if combinations.size == 0:
            break
        energies = fun(combinations)
        minimum = np.argmin(energies)
        if energies[minimum] < energy:
            energy = energies[minimum]
            result = combinations[:, minimum]
    y = int("".join(str(int(e)) for e in result[:H.ny]) + "1", base=2)
    x = int("".join(str(int(e)) for e in result[H.ny:H.ny+H.nx]) + "1", base=2)
    return energy, x, y

def main():
    for N in range(33, 90, 2):
        H = DirectHamiltonianSympy(N).with_parameters(letters=False, show_products=True)
        print(H)
        # exp0 = str(sympy.expand(H.from_simpy_formula()))
        # exp1 = str(sympy.expand(H.from_alternate_formula()))
        # sym1, res1, x1, y1 = minimize(exp1, H)
        # assert(exp0 == exp1)

        # exp2 = str(sympy.expand(H.from_ishikawa_formula()))
        # print(exp2)
        # sym2, res2, x2, y2 = minimize(exp2, H)
        # print(f"[ISHI] x={x2:3d} y={y2:3d} E={res2[0]}")
        # assert(res1[0] == res2[0])
        # assert(x1 == x2)
        # assert(y1 == y2)

        H2 = DirectHamiltonian(N).with_parameters(letters=False, show_products=True)
        exp2 = str(sympy.expand(H2.quadratized_by_dattani_chau()))
        # print(exp2)

        sym2, res2, x2, y2 = minimize(exp2, H2)
        print(f"[DATT_OLD] x={x2:3d} y={y2:3d} E={res2[0]}")
        # assert(res1[0] == res2[0])
        # assert(x1 == x2)
        # assert(y1 == y2)

        exp3 = H.from_dattani_formula()
        exp3 = str(sympy.expand(exp3))
        # print(exp3)
        sym3, res3, x3, y3 = minimize(exp3, H)
        print(f"[DATT_NEW] x={x3:3d} y={y3:3d} E={res3[0]}")
        assert(exp2 == exp3)
        # print(sym3)

        # assert(res1[0] == res3[0])
        # assert(x1 == x3)
        # assert(y1 == y3)

def sanitize(string):
    string = string.replace("*", "")
    string = string.replace("x0", "x")
    string = string.replace("y0", "y")
    string = string.replace("b0", "b")
    return string

def test():
    for N in range(25, 80, 1):
        H = DirectHamiltonianSympy(N).with_parameters(letters=False, show_products=True)
        print(H)
        ishi = str(sympy.expand(H.from_ishikawa_formula()))
        if N <= 50:
            energy_ishi, x_ishi, y_ishi = minimize(ishi, H)
            print(f"[ISHI] x={x_ishi:3d}; y={y_ishi:3d}; E={energy_ishi:3d}")
        else:
            print("[ISHI] Can't do it!")
        # print(sanitize(ishi), end="\n")

        datt = str(sympy.expand(H.from_dattani_formula()))
        if N <= 80:
            energy_datt, x_datt, y_datt = minimize(datt, H)
            print(f"[DATT] x={x_datt:3d}; y={y_datt:3d}; E={energy_datt:3d}")
        else:
            print("[DATT] Can't do it!")
        # print(sanitize(datt))

if __name__ == "__main__":
    # main()
    test()
