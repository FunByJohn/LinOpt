from linear_program import *
from matrix import *

def primal_simplex(LP, initial_basis):
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

    print(tableau)


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