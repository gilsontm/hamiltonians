import sympy
from direct_hamiltonian import DirectHamiltonian


def minimize(expression, H):
    exp = sympy.sympify(expression)
    symbols = list(exp.free_symbols)
    symbols.sort(key=lambda s: s.name, reverse=True)
    fun = sympy.lambdify(symbols, exp, "numpy")
    attrs = [bin(x)[2:].rjust(len(symbols), '0') for x in range(2**(len(symbols)))]
    results = [(fun(*map(int, attr)), attr) for attr in attrs]
    results.sort()
    y = int(results[0][1][:H.ny] + '1', base=2)
    x = int(results[0][1][H.ny:H.ny+H.nx] + '1', base=2)
    return symbols, results, x, y

def main():
    Ns = [25, 33, 35, 39, 49]
    for N in range(25, 50, 2):
        H = DirectHamiltonian(N).with_parameters(letters=False, show_products=True)
        print(H)
        exp0 = H.from_function()
        exp0 = sympy.expand(exp0)
        exp0 = str(exp0).replace("**2", "")
        exp0 = str(sympy.simplify(exp0))
        symbols, res, x, y = minimize(exp0, H)
        print(f"[Orig] x = {x:3d}, y = {y:3d}, x*y = {x*y:3d}")
        # print((res[0][0], res[0][1][:len(symbols)]))

        exp2 = H.quadratized_by_ishikawa()
        exp2 = str(sympy.simplify(exp2))
        sym, res, x, y = minimize(exp2, H)
        print(f"[Ishi] x = {x:3d}, y = {y:3d}, x*y = {x*y:3d}")
        # print((res[0][0], res[0][1][:len(symbols)]))

        exp3 = H.quadratized_by_dattani_chau()
        exp3 = str(sympy.simplify(exp3))
        sym, res, x, y = minimize(exp3, H)
        print(f"[Datt] x = {x:3d}, y = {y:3d}, x*y = {x*y:3d}")
        # print((res[0][0], res[0][1][:len(symbols)]))

if __name__ == "__main__":
    main()
