import ipdb
import __init__
import numpy as np
from utils.utils import *
import matplotlib.pyplot as plt
from number_of_bits import number_of_bits_for_x, number_of_bits_for_y


# class Term:
#     def __init__(self, multiplier, string):
#         self.multiplier = multiplier
#         self.string = string

#     def combine(self, other):
#         return Term(self.multiplier * other.multiplier, self.string + other.string)

#     def __str__(self):
#         return f"{self.multiplier}*{self.string}"

def expression_for_sat(clauses):
    s = ""
    for clause in clauses:
        for item in clause:
            char = chr(96 + abs(item))
            char = "x" if char == "e" else char
            if item > 0:
                s += f"(1-{char})"
            else:
                s += f"{char}"
        s += " + "
    return s[:-3]

def expressions(N, nx, ny):
    s = f"{N:3d},{nx:3d},{ny:3d} -> ({N - 1:2d} - "
    for i in range(1, nx + 1):
        for j in range(1, ny + 1):
            char1 = chr(96+i)
            char1 = "x" if char1 == "e" else char1
            char2 = chr(96 + nx + j)
            char2 = "x" if char2 == "e" else char2
            s += f"{2**(i+j):2d}{char1}{char2} - "
    for i in range(1, nx + 1):
        char2 = chr(96 + i)
        char2 = "x" if char2 == "e" else char2
        s += f"{2**i:2d}{char2} - "
    for j in range(1, ny + 1):
        char2 = chr(96 + nx + j)
        char2 = "x" if char2 == "e" else char2
        s += f"{2**j:2d}{char2} - "
    s = s[:-3] + ")^2"
    return s


def main():
    # numbers = range(9, 130, 1)
    numbers = [33]
    xs = map(lambda N: int(number_of_bits_for_x(N)), numbers)
    ys = map(lambda N: int(number_of_bits_for_y(N)), numbers)
    tuples = list(zip(numbers, xs, ys))
    for t in tuples:
        print(expressions(*t))

    # clauses = [(1, 2, 3, 4), (-2, -3, 4, 5), (1, 2, 3, 5), (-1, 2, 4, 5)]

    # clauses = [(-1,-1,-2,-2,),
    #             (-1,-1,-2,-3,),
    #             (-1,-1,-3,-3,),
    #             (-1,-1,-2,),
    #             (-1,-1,-3,),
    #             (-1,-2,-2,),
    #             (-1,-2,-3,),
    #             (-1,-3,-3,),
    #             (-1,-1,),
    #             (-1,-2,),
    #             (-1,-3,),
    #             (-2,-2,),
    #             (-2,-3,),
    #             (-3,-3,),
    #             (-1,),
    #             (-2,),
    #             (-3,),]

    # print(expression_for_sat(clauses))

if __name__ == "__main__":
    main()