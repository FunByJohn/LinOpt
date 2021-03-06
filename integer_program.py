import random

from linear_program import *
from simplex import *

class IntegerProgram(LinearProgram):

    def __init__(self):
        super().__init__()

    def is_integer_feasible(self, vector):
        if not super().is_feasible(vector):
            return False

        # check if integer
        for value in vector:
            if not value.is_integer():
                return False

        return True

# Bug: Sometimes it breaks the dual algorithm, just run it again and it might work
def randomized_gomory(IP, printTeX = True):
    iterations = 0

    # Algorithm on p. 482 in the book
    relaxedLP = IP.get_standard_form()
    initial_basis = find_initial_basic_solution(relaxedLP, printTeX = printTeX)

    if initial_basis == None:
        print("Error in Randomized Gomory: Couldn't find initial basic feasible solution to relaxed LP!")
        return None

    sol, basis, tableau = simplex(relaxedLP, initial_basis, SimplexAlg.PRIMAL, printTeX = printTeX)
    orig_sol = sol[0:len(IP.variable_names)]

    while not IP.is_integer_feasible(orig_sol):
        non_basic = []
        for i in range(len(relaxedLP.variable_names)):
            if i not in basis:
                non_basic.append(i)

        candidates = []
        for i in range(1, tableau.m):
            value = tableau.at(i, 0)
            if not value.is_integer():
                candidates.append((i - 1, Q_floor(value)))

        (k, value) = random.choice(candidates)

        lhs = []
        lhs.append((Q(1), relaxedLP.variable_names[k]))
        
        for j in non_basic:
            coef = Q_floor(tableau.at(k + 1, j + 1))
            if coef.is_nonzero():
                lhs.append((coef, relaxedLP.variable_names[j]))

        if printTeX:
            print(f'Added the Gomory cut ${TeXify_linear_combination(lhs)} \\leq {value.to_TeX()}$')

        relaxedLP.add_constraint([lhs, ComparisonOperator.LESS_EQUAL, value])
        relaxedLP = relaxedLP.get_standard_form()

        basis.append(len(relaxedLP.variable_names) - 1)

        sol, basis, tableau = simplex(relaxedLP, basis, SimplexAlg.DUAL, printTeX = printTeX)
        orig_sol = sol[0:len(IP.variable_names)]

        iterations += 1

    if printTeX:
        print(f'Randomized Gomory used {iterations} iterations this time!')

    return orig_sol
