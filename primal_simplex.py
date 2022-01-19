from linear_program import *
from matrix import *

import math

def primal_simplex(LP, initial_basis, printTeX = True):
    assert LP.standard_form == True

    current_basis = list(map(lambda x: LP.variable_names.index(x), initial_basis.split(' ')))

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

    def get_tableau_TeX(tableau):
        '''
        Example:

        \bgroup
        \def\arraystretch{1.5}
        \begin{table}[H]
            \begin{tabular}{r|r|rrrrrr|}
            \cline{2-8}
                    &       & $x_1$  & $x_2$  & $x_3$  & $x_4$  & $x_5$  & $x_6$  \\ \cline{2-8} 
                    & $-10$ & $0$    & $1$    & $-4$   & $-2$   & $0$    & $0$    \\ \cline{2-8} 
            $x_1 =$ & $5$   & $1$    & $0$    & $3$    & $1$    & $0$    & $0$    \\
            $x_5 =$ & $3$   & $0$    & $2$    & $1$    & $0$    & $1$    & $0$    \\
            $x_6 =$ & $8$   & $0$    & $-1$   & $5$*   & $2$    & $0$    & $1$    \\ \cline{2-8}
            \end{tabular}
        \end{table}
        \egroup
        '''

        TeX = []
        TeX.append(f'\\bgroup')
        TeX.append(f'\\def\\arraystretch{{1.5}}')
        TeX.append(f'\\begin{{table}}[H]')

        rs = f'r|r|' + ('r' * A.n) + '|'
        cline = f'\\cline{{2-{str(A.n + 2)}}}'
        
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

    print(get_tableau_TeX(tableau))

    # Now we follow the algorithm as described on p. 100
    # Step 1 has been performed in the code above
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

    solution = [Q(0) for i in range(len(LP.variable_names))]

    for i in range(len(current_basis)):
        solution[current_basis[i]] = tableau.at(i + 1, 0)

    return solution

LP = LinearProgram()
LP.set_objective('  min  2 x_1   + x_2 + 2 x_3')
LP.add_constraint('        x_1         + 3 x_3  <= 5 ')
LP.add_constraint('              2 x_2 +   x_3  <= 3 ')
LP.add_constraint('      2 x_1 +   x_2 +   x_3  >= 2 ')
#LP.add_constraint('        x_1                  >= 0 ')
#LP.add_constraint('                x_2          >= 0 ') # TODO: How to handle x_1, x_2, x_3 >= 0 constraint?
#LP.add_constraint('                        x_3  >= 0 ')

LP_standard = LP.get_standard_form()
primal_simplex(LP_standard, 'x_1 x_5 x_6')