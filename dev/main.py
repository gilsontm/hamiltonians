import sympy
from direct_hamiltonian import DirectHamiltonian

if __name__ == "__main__":
    # for N in range(11, 30, 4):
    N = 33
    H = DirectHamiltonian(N).with_parameters(letters=False, show_products=True)
    print(str(H))
    exp0 = H.from_expansion()
    print(f"H = {exp0}")
    if H.show_products:
        exp0 = str(sympy.simplify(exp0))
        print(f"H = {exp0}\n")

    print("Ishikawa:")
    exp2 = H.quadratized_by_ishikawa()
    print(f"H = {exp2}")
    if H.show_products:
        exp2 = str(sympy.simplify(exp2))
        print(f"H = {exp2}\n")

    print("Dattani-Chau:")
    exp1 = H.quadratized_by_dattani_chau()
    print(f"H = {exp1}")
    if H.show_products:
        exp1 = str(sympy.simplify(exp1))
        print(f"H = {exp1}\n")

        # exp0 = H.from_function()
        # exp0 = sympy.expand(exp0)
        # exp0 = str(exp0).replace("**2", "")
        # exp0 = str(sympy.simplify(exp0))
        # exp1 = str(sympy.simplify(H.from_expansion()))
        # exp2 = str(sympy.simplify(H.from_formula()))
        # assert(exp0 == exp1 and exp0 == exp2)
        # print(exp2)