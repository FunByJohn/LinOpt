from enum import Enum, auto
from linear_program import *
from matrix import *

import math

def get_tableau_TeX(LP, current_basis, tableau):
    TeX = []
    TeX.append(f'\\bgroup')
    TeX.append(f'\\def\\arraystretch{{1.5}}')
    TeX.append(f'\\begin{{table}}[H]')

    rs = f'r|r|' + ('r' * (tableau.n - 1)) + '|'
    cline = f'\\cline{{2-{str((tableau.n - 1) + 2)}}}'
    
    TeX.append(f'  \\begin{{tabular}}{{{rs}}}')
    TeX.append(f'  {cline}')
    
    xs = ' & '.join([f'${x}$' for x in LP.variable_names])
    TeX.append(f'    &   & {xs} \\\\ {cline}')

    cs = ' & '.join([f'${tableau.at(0, i).to_TeX()}$' for i in range(tableau.n)])
    TeX.append(f'    & {cs} \\\\ {cline}')

    for i in range(1, tableau.m):
        the_cline = ''  if i < tableau.m - 1 else cline
        vs = ' & '.join([f'${tableau.at(i, j).to_TeX()}$' for j in range(tableau.n)])
        TeX.append(f'  ${LP.variable_names[current_basis[i-1]]} = $ & {vs} \\\\ {the_cline}')
    
    TeX.append(f'  \\end{{tabular}}')
    TeX.append(f'\\end{{table}}')
    TeX.append(f'\\egroup')

    return '\n'.join(TeX)

# Needs to rows of A to be lineaarly independent to work, but if A with linearly dependent rows
# is supplied, the algorithm will detect this and let us now
def find_initial_basic_solution(LP, printTeX = True):
    # We follow the algorithm on p. 116-117
    LP_aux = LP.get_standard_form()

    # Step 1
    for i, (lhs, comparison, rhs) in enumerate(LP_aux.constraints):
        if rhs.is_negative():
            new_lhs = []
            new_rhs = rhs.add_inv()

            for (coef, variable) in lhs:
                new_lhs.append((coef.add_inv(), variable))

            LP_aux.constraints[i] = (new_lhs, comparison, new_rhs)

    # Step 2
    LP_aux.objective = []

    for i, (lhs, comparison, rhs) in enumerate(LP_aux.constraints):
        new_var = LP_aux.get_new_variable_name()
        lhs.append((Q(1), new_var))
        LP_aux.objective.append((Q(1), new_var))

    if printTeX:
        print("First we set up the auxillary linear program")
        print(LP_aux.get_TeX())

    # Step 3
    initial_basis = [i for i in range(len(LP_aux.variable_names) - len(LP_aux.constraints), len(LP_aux.variable_names))]
    solution, new_basis, tableau = simplex(LP_aux, initial_basis, SimplexAlg.PRIMAL, printTeX = printTeX)

    if tableau.at(0, 0).is_positive():
        print("Couldn't find initial basic solution as the problem is infeasible!")
        return None

    if printTeX:
        print(f"Optimal solution is {tableau.at(0, 0)}, so the problem is feasible, and we now find the initial basis")

    # Step 4
    # Optimal cost is now zero, since x_i >= 0 for all i, so the auxillary objective function can only take nonnegative values
    def get_artificial_var_in_basis():
        for i in new_basis:
            if i in initial_basis:
                return i
        return None

    if get_artificial_var_in_basis() == None:
        if printTeX:
            print("None of the new variables are in the basis, so we have found an initial basic feasible solution!")
        return new_basis

    # Step 5
    while True:
        l = get_artificial_var_in_basis()

        if l == None:
            break

        l = new_basis.index(l)

        # Detect if a row is redundant
        is_redundant = True
        
        j = -1
        for i in range(1, len(LP.variable_names) + 1):
            if not tableau.at(l + 1, i).is_zero():
                is_redundant = False
                j = i - 1
                break

        if is_redundant:
            print("Original LP has a redundant constraint, that is the rows of the original A are not linearly independent!")
            print("Find RREF of A^T (for the original A) in order to see which constraint is redundant, and drop it before calling this function.")
            assert False

        # Change basis to eliminate artificial variable
        new_basis[l] = j
        tableau = tableau.scale_row(l+1, tableau.at(l + 1, j + 1).mul_inv())

        for i in range(0, tableau.m):
            if i != l + 1:
                tableau = tableau.add_rows(i, l + 1, tableau.at(i, j + 1).add_inv())

        if printTeX:
            print(get_tableau_TeX(LP_aux, new_basis, tableau))

    if printTeX:
        print("None of the new variables are in the basis, so we have found an initial basic feasible solution!")
    
    return new_basis

class SimplexAlg(Enum):
    PRIMAL = auto()
    DUAL = auto()

def simplex(LP, initial_basis, alg, printTeX = True):
    assert LP.standard_form == True

    if isinstance(initial_basis, str):
        current_basis = list(map(lambda x: LP.variable_names.index(x), initial_basis.split(' ')))
    else:
        current_basis = initial_basis.copy()

    def make_A():
        rows = []
        
        for (lhs, comparison, rhs) in LP.constraints:
            row = [Q(0) for i in range(len(LP.variable_names))]

            for (coef, variable) in lhs:
                row[LP.variable_names.index(variable)] = coef

            rows.append(row)

        return Matrix(*rows)

    def pick_columns_of_current_basis(A):
        columns = []

        for i in current_basis:
            columns.append(A.submatrix(i, 0, 1, A.m))

        B = columns[0]

        for i in range(1, len(columns)):
            B = Matrix_combine_horizontal(B, columns[i])

        return B

    def make_b():
        coefs = [Q(0) for i in LP.constraints]

        for i in range(len(LP.constraints)):
            (lhs, comparison, rhs) = LP.constraints[i]
            coefs[i] = rhs

        return Matrix(*[coefs]).transpose()

    def make_c():
        coefs = [Q(0) for i in LP.variable_names]

        for (coef, variable) in LP.objective:
            coefs[LP.variable_names.index(variable)] = coef

        return Matrix(*[coefs]).transpose()

    A = make_A()
    B = pick_columns_of_current_basis(A)
    b = make_b()
    c = make_c()
    B_inv = B.inverse()
    c_B = pick_columns_of_current_basis(c.transpose()).transpose()

    top_left = Matrix_mul(Matrix(*[[Q(-1)]]), Matrix_mul(Matrix_mul(c_B.transpose(), B_inv), b))
    top_right = Matrix_sub(c.transpose(), Matrix_mul(Matrix_mul(c_B.transpose(), B_inv), A))
    bot_left = Matrix_mul(B_inv, b)
    bot_right = Matrix_mul(B_inv, A)
    tableau = Matrix_combine_vertical(Matrix_combine_horizontal(top_left, top_right), Matrix_combine_horizontal(bot_left, bot_right))

    def current_solution():
        solution = [Q(0) for i in range(len(LP.variable_names))]

        for i in range(len(current_basis)):
            solution[current_basis[i]] = tableau.at(i + 1, 0)

        return solution

    if printTeX:
        print(get_tableau_TeX(LP, current_basis, tableau))

    #print('===============')
    #print(Matrix_rref(A.transpose()).transpose())
    #print('===============')

    if alg == SimplexAlg.PRIMAL:
        # Now we follow the algorithm as described on p. 100
        
        # Step 1
        if not LP.is_feasible(current_solution()):
            print('WARNING! Initial solution is not feasible')
            assert False

        while (True):
            # Step 2
            is_optimal = True
            j = 0
            
            for i in range(1, tableau.n):
                if tableau.at(0, i).is_negative():
                    j = i - 1
                    is_optimal = False
                    break

            if is_optimal:
                break

            # Step 3
            u = tableau.submatrix(j + 1, 1, 1, A.m)

            has_positive = False
            for i in range(u.m):
                if u.at(i, 0).is_positive():
                    has_positive = True
                    break

            if not has_positive:
                print("Optimal solution is -infinity")
                break

            # Step 4
            l = -1
            smallest_ratio = math.inf

            for i in range(u.m):
                if u.at(i, 0).is_positive():
                    ratio = Q_div(tableau.at(i + 1, 0), u.at(i, 0))

                    if l == -1 or Q_leq(ratio, smallest_ratio):
                        l = i
                        smallest_ratio = ratio

            # print(f'{LP.variable_names[current_basis[l]]} leaves the basis')
            # print(f'{LP.variable_names[j]} enters the basis')

            current_basis[l] = j

            # Step 5
            tableau = tableau.scale_row(l+1, tableau.at(l + 1, j + 1).mul_inv())

            for i in range(0, tableau.m):
                if i != l + 1:
                    tableau = tableau.add_rows(i, l + 1, tableau.at(i, j + 1).add_inv())

            if printTeX:
                print(get_tableau_TeX(LP, current_basis, tableau))
    
    elif alg == SimplexAlg.DUAL:
        # Now we follow the algorithm as described on p. 159

        # Step 1
        for i in range(1, tableau.n):
            if tableau.at(0, i).is_negative():
                print('WARNING! Initial solution cannot be used by dual simplex, a reduced cost is negative!')
                assert False

        while (True):
            l = -1

            # Step 2
            is_optimal = True
            for i in range(1, tableau.m):
                if tableau.at(i, 0).is_negative():
                    is_optimal = False
                    l = i - 1

            if is_optimal:
                break

            # Step 3
            is_optimal = True

            for i in range(1, tableau.n):
                if tableau.at(l + 1, i).is_negative():
                    is_optimal = False
                    break

            if is_optimal:
                print('Optimal solution is infinity')
                break

            # Step 4
            j = math.inf
            smallest_ratio = math.inf

            for i in range(1, tableau.n):
                if tableau.at(l + 1, i).is_negative():
                    ratio = Q_div(tableau.at(0, i), Q_abs(tableau.at(l + 1, i)))

                    if j == math.inf or Q_leq(ratio, smallest_ratio):
                        j = i - 1
                        smallest_ratio = ratio

            #print(f'{LP.variable_names[current_basis[l]]} leaves the basis')
            #print(f'{LP.variable_names[j]} enters the basis')

            current_basis[l] = j

            # Step 5
            tableau = tableau.scale_row(l + 1, tableau.at(l + 1, j + 1).mul_inv())

            for i in range(0, tableau.m):
                if i != l + 1:
                    tableau = tableau.add_rows(i, l + 1, tableau.at(i, j + 1).add_inv())

            if printTeX:
                print(get_tableau_TeX(LP, current_basis, tableau))

    solution = current_solution()

    if printTeX:
        print('Optimal solution is')
        print('\\[')

        variables = ', '.join(LP.variable_names)
        values = ', '.join([value.to_TeX() for value in solution])
        print(f'  ({variables}) = ({values})')
        print('\\]')

    if not LP.is_feasible(solution):
        print("WARNING! Found solution is not feasible")
        assert False

    return solution, current_basis, tableau