import random

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

def randomized_gomory(IP, printTeX = True):
    iterations = 0

    # Algorithm on p. 482 in the book
    relaxedLP = IP.get_standard_form()

    # Solve relaxed LP form of IP   (TODO: Make this more general)
    sol, basis, tableau = simplex(relaxedLP, [0, 1], SimplexAlg.PRIMAL, printTeX = printTeX)
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

#IP = IntegerProgram()
#IP.set_objective('min x_1 - 2x_2')
#IP.add_constraint('-4x_1 + 6x_2 <= 9')
#IP.add_constraint('x_1 + x_2 <= 4')

#print(randomized_gomory(IP, True))

IP = IntegerProgram()
IP.set_objective('min x_2 + x_3')
IP.add_constraint('x_1 + x_2 + x_3 >= 4')
IP.add_constraint('4x_1 + 2x_2 + x_3 <= 12')

print(randomized_gomory(IP, True))