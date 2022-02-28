import sympy
import itertools
import numpy as np
from labels import *
from direct_hamiltonian_sympy import DirectHamiltonianSympy


class Tester:

    def __init__(self, verbose):
        self.verbose = verbose

    @staticmethod
    def __sympify(expression):
        sympified = sympy.sympify(expression)
        sympified = sympy.expand(sympified)
        symbols = list(sympified.free_symbols)
        symbols.sort(key=lambda symbol: symbol.name, reverse=True)
        function = sympy.lambdify((symbols,), sympified, "numpy")
        parameters = len(symbols)
        return function, parameters

    @staticmethod
    def __minimize(expression, hamiltonian):
        WINDOW = 256
        energy = 1e10
        ny = hamiltonian.ny
        nx = hamiltonian.nx

        function, parameters = Tester.__sympify(expression)
        iterator = itertools.product([0, 1], repeat=parameters)

        while True:
            combinations = np.bool_(list(itertools.islice(iterator, WINDOW))).transpose()
            if combinations.size == 0:
                break
            energies = function(combinations)
            minimum = np.argmin(energies)
            if energies[minimum] < energy:
                energy = energies[minimum]
                result = combinations[:, minimum]
        y = int("".join(str(int(e)) for e in result[:ny]) + "1", base=2)
        x = int("".join(str(int(e)) for e in result[ny:ny+nx]) + "1", base=2)
        return energy, x, y

    @staticmethod
    def __minimize_partially(expression, hamiltonian, prefix):
        WINDOW = 256
        energy = 1e10
        ny = hamiltonian.ny
        nx = hamiltonian.nx
        prefix = np.bool_(prefix)
        function, parameters = Tester.__sympify(expression)
        iterator = itertools.product([0, 1], repeat=(parameters - nx - ny))

        while True:
            combinations = np.bool_(list(itertools.islice(iterator, WINDOW)))
            if combinations.size == 0:
                break
            tiled = np.tile(prefix, (combinations.shape[0], 1))
            combinations = np.concatenate((tiled, combinations), axis=1).transpose()
            energies = function(combinations)
            minimum = np.argmin(energies)
            if energies[minimum] < energy:
                energy = energies[minimum]
                result = combinations[:, minimum]
        y = int("".join(str(int(e)) for e in result[:ny]) + "1", base=2)
        x = int("".join(str(int(e)) for e in result[ny:ny+nx]) + "1", base=2)
        return energy, x, y

    def test_original_form(self):
        FROM = 25; TO = 1024
        print(TEST_0.format("Original vs. Alternate Form"))
        print(FROM_0_THROUGH_1.format(FROM, TO))
        for N in range(FROM, TO, 2):
            H = DirectHamiltonianSympy(N)
            original = H.from_simpy_formula()
            original = sympy.sympify(original)
            alternate = H.from_alternate_formula()
            alternate = sympy.sympify(alternate)
            assert original == alternate, ERROR_AT_N_0.format(N)
            if self.verbose:
                print(PASSED_FOR_N_0.format(N))
        print(PASSED_ALL)

    def test_ishikawa_minimization(self):
        FROM = 25; TO = 50
        print(TEST_0.format("Ishikawa Minimization"))
        print(FROM_0_THROUGH_1.format(FROM, TO))

        for N in range(FROM, TO, 2):
            H = DirectHamiltonianSympy(N)
            original = H.from_simpy_formula()
            ishikawa = H.from_ishikawa_formula()

            energy_orig, x_orig, y_orig = Tester.__minimize(original, H)
            energy_ishi, x_ishi, y_ishi = Tester.__minimize(ishikawa, H)

            assert energy_orig == energy_ishi, ERROR_AT_N_0_EXPECTED_1_BUT_WAS_2.format(N, energy_orig, energy_ishi)
            assert x_orig == x_ishi,           ERROR_AT_N_0_EXPECTED_1_BUT_WAS_2.format(N, x_orig, x_ishi)
            assert y_orig == y_ishi,           ERROR_AT_N_0_EXPECTED_1_BUT_WAS_2.format(N, y_orig, y_ishi)

            if self.verbose:
                print(AT_N_0_ENERGY_1_X_2_Y_3.format(N, energy_ishi, x_ishi, y_ishi))
        print(PASSED_ALL)

    def test_dattani_minimization(self):
        FROM = 25; TO = 80
        print(TEST_0.format("Dattani-Chau Minimization"))
        print(FROM_0_THROUGH_1.format(FROM, TO))

        for N in range(FROM, TO, 2):
            H = DirectHamiltonianSympy(N)
            original = H.from_simpy_formula()
            dattani = H.from_dattani_formula()

            energy_orig, x_orig, y_orig = Tester.__minimize(original, H)
            energy_datt, x_datt, y_datt = Tester.__minimize(dattani, H)

            assert energy_orig == energy_datt, ERROR_AT_N_0_EXPECTED_1_BUT_WAS_2.format(N, energy_orig, energy_datt)
            assert x_orig == x_datt,           ERROR_AT_N_0_EXPECTED_1_BUT_WAS_2.format(N, x_orig, x_datt)
            assert y_orig == y_datt,           ERROR_AT_N_0_EXPECTED_1_BUT_WAS_2.format(N, y_orig, y_datt)

            if self.verbose:
                print(AT_N_0_ENERGY_1_X_2_Y_3.format(N, energy_datt, x_datt, y_datt))
        print(PASSED_ALL)

    def test_factored_values_evaluation(self):
        FROM = 25; TO = 1024
        print(TEST_0.format("Factored Values Evaluation"))
        print(FROM_0_THROUGH_1.format(FROM, TO))

        for N in range(FROM, TO, 2):
            factors = sympy.ntheory.factorint(N)
            if len(factors.values()) != 2 or sum(factors.values()) != 2:
                continue
            H = DirectHamiltonianSympy(N)
            original = H.from_simpy_formula()
            alternate = H.from_simpy_formula()
            function_orig, _ = Tester.__sympify(original)
            function_alte, _ = Tester.__sympify(alternate)

            x, y = sorted(list(factors.keys()))
            x_bin = tuple(map(int, bin(x)[2:][:-1].rjust(H.nx, "0")))
            y_bin = tuple(map(int, bin(y)[2:][:-1].rjust(H.ny, "0")))
            combination = y_bin + x_bin

            energy_orig = function_orig(combination)
            energy_alte = function_alte(combination)

            assert energy_orig == 0, ERROR_AT_N_0_EXPECTED_1_BUT_WAS_2.format(N, 0, energy_orig)
            assert energy_alte == 0, ERROR_AT_N_0_EXPECTED_1_BUT_WAS_2.format(N, 0, energy_alte)
            if self.verbose:
                print(AT_N_0_ENERGY_1_X_2_Y_3.format(N, 0, x, y))
        print(PASSED_ALL)

    def test_factored_values_evaluation_with_ishikawa_minimization(self):
        FROM = 25; TO = 80
        print(TEST_0.format("Factored Values Evaluation with Ishikawa Minimization"))
        print(FROM_0_THROUGH_1.format(FROM, TO))

        for N in range(FROM, TO, 2):
            factors = sympy.ntheory.factorint(N)
            if len(factors.values()) != 2 or sum(factors.values()) != 2:
                continue
            H = DirectHamiltonianSympy(N)
            x, y = sorted(list(factors.keys()))
            x_bin = tuple(map(int, bin(x)[2:][:-1].rjust(H.nx, "0")))
            y_bin = tuple(map(int, bin(y)[2:][:-1].rjust(H.ny, "0")))
            prefix = y_bin + x_bin

            ishikawa = H.from_ishikawa_formula()
            energy_ishi, _, _ = Tester.__minimize_partially(ishikawa, H, prefix)

            assert energy_ishi == 0, ERROR_AT_N_0_EXPECTED_1_BUT_WAS_2.format(N, 0, energy_ishi)
            if self.verbose:
                print(AT_N_0_ENERGY_1_X_2_Y_3.format(N, energy_ishi, x, y))
        print(PASSED_ALL)

    def test_factored_values_evaluation_with_dattani_minimization(self):
        FROM = 25; TO = 95
        print(TEST_0.format("Factored Values Evaluation with Dattani Minimization"))
        print(FROM_0_THROUGH_1.format(FROM, TO))

        for N in range(FROM, TO, 2):
            factors = sympy.ntheory.factorint(N)
            if len(factors.values()) != 2 or sum(factors.values()) != 2:
                continue
            H = DirectHamiltonianSympy(N)
            x, y = sorted(list(factors.keys()))
            x_bin = tuple(map(int, bin(x)[2:][:-1].rjust(H.nx, "0")))
            y_bin = tuple(map(int, bin(y)[2:][:-1].rjust(H.ny, "0")))
            prefix = y_bin + x_bin

            dattani = H.from_dattani_formula()
            energy_datt, _, _ = Tester.__minimize_partially(dattani, H, prefix)

            assert energy_datt == 0, ERROR_AT_N_0_EXPECTED_1_BUT_WAS_2.format(N, 0, energy_datt)
            if self.verbose:
                print(AT_N_0_ENERGY_1_X_2_Y_3.format(N, 0, x, y))
        print(PASSED_ALL)

    @staticmethod
    def run(verbose=False):
        tester = Tester(verbose=verbose)
        tester.test_original_form()
        tester.test_ishikawa_minimization()
        tester.test_dattani_minimization()
        tester.test_factored_values_evaluation()
        tester.test_factored_values_evaluation_with_ishikawa_minimization()
        tester.test_factored_values_evaluation_with_dattani_minimization()


if __name__ == "__main__":
    Tester.run(True)
