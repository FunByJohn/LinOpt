from enum import Enum, auto
from rational import *
from matrix import *
from simplex import *

import re
import copy

def TeXify_linear_combination(lin_comb):
    string = ''

    for i in range(len(lin_comb)):
        (coef, variable) = lin_comb[i]

        if coef.is_one():
            coef_str = ''
        elif Q_mul(coef, Q(-1)).is_one():
            coef_str = '-'
        else:
            coef_str = coef.to_TeX()

        if i > 0:
            if not coef.is_negative():
                string += '+'

        string += coef_str + variable

    return string

class MinOrMax(Enum):
    MIN = auto()
    MAX = auto()

class ComparisonOperator(Enum):
    EQUAL         = auto()
    LESS_EQUAL    = auto()
    GREATER_EQUAL = auto()

class LinearProgram:

    def __init__(self):
        self.min_or_max = None
        self.objective = []
        self.variable_names = []
        self.constraints = []
        self.standard_form = True
        self.is_dual = False

    def set_objective(self, string):
        (min_max, lin_comb) = string.lstrip().rstrip().split(' ', maxsplit = 1)

        self.min_or_max = MinOrMax.MIN
        self.objective = self.parse_linear_combination(lin_comb)

        if min_max == 'max':
            # We flip signs in objective so it becomes min instead, since set_objective will only be used on primal problems
            new_objective = []

            for (coef, variable) in self.objective:
                new_objective.append((coef.add_inv(), variable))

            self.objective = new_objective
        elif min_max != 'min':
            print('Error when initializing LinearProgram: "min" or "max" not supplied in objective function')
            assert False

    def add_constraint(self, string):
        if isinstance(string, str):
            comparisonStr = re.search(r'([<>]{0,1}=)', string).group(1)
            (lhs, rhs) = re.split(r'\s*[<>]{0,1}=\s*', string)
            lhs = self.parse_linear_combination(lhs)
            rhs = Q(rhs)

            if comparisonStr == '=':
                comparison = ComparisonOperator.EQUAL
            else:
                self.standard_form = False
                
                if comparisonStr == '<=':
                    comparison = ComparisonOperator.LESS_EQUAL
                elif comparisonStr == '>=':
                    comparison = ComparisonOperator.GREATER_EQUAL
        else:
            (lhs, comparison, rhs) = string

        for (coef, variable) in lhs:
            self.ensure_variable_exists(variable)

        self.constraints.append([lhs, comparison, rhs])

    def ensure_variable_exists(self, string):
        if not string in self.variable_names:
            self.variable_names.append(string)
            self.variable_names = sorted(self.variable_names, key = lambda name: int(name.split('_')[1]))

    def parse_linear_combination(self, string):
        result = []
        parts = re.findall(r'[-\+]{0,1}\s*\d*\s*[a-z]_\d+', string)

        for part in parts:
            match = re.search(r'([-\+]{0,1}\s*\d*)\s*([a-z]_\d+)', part)
            coef = match.group(1).lstrip().rstrip()
            name = match.group(2)
            sgn = 1 if re.match(r'-', coef) == None else -1

            if coef == '-':
                coef = Q(-1)
            else:
                coef = re.sub(r'[-\+]', '', coef)

                if coef == '':
                    coef = Q(1)
                else:
                    coef = Q_mul(Q(sgn), Q(coef))

            result.append((coef, name))

        return result

    def is_feasible(self, vector):
        assert len(vector) == len(self.variable_names)

        if self.is_dual == False:
            for value in vector:
                if value.is_negative():
                    return False

        for (lhs, comparison, rhs) in self.constraints:
            lhs_value = Q(0)

            for (coef, variable) in lhs:
                lhs_value = Q_add(lhs_value, Q_mul(coef, vector[self.variable_names.index(variable)]))

            if comparison == ComparisonOperator.EQUAL:
                if not Q_equals(lhs_value, rhs):
                    return False
            elif comparison == ComparisonOperator.LESS_EQUAL:
                if not Q_leq(lhs_value, rhs):
                    return False
            elif comparison == ComparisonOperator.GREATER_EQUAL:
                if not Q_geq(lhs_value, rhs):
                    return False

        return True

    def get_new_variable_name(self):
        last_variable_index = int(self.variable_names[-1].split('_')[1])
        new_variable = f'x_{last_variable_index + 1}'
        self.ensure_variable_exists(new_variable)
        return new_variable

    # slack and surplus
    def add_artificial_variables(self):
        self.standard_form = True

        for i in range(len(self.constraints)):
            (lhs, comparison, rhs) = self.constraints[i]

            if comparison == ComparisonOperator.LESS_EQUAL:
                # add slack variable
                lhs.append((Q(1), self.get_new_variable_name()))
            elif comparison == ComparisonOperator.GREATER_EQUAL:
                # add surplus variable
                lhs.append((Q(-1), self.get_new_variable_name()))

            comparison = ComparisonOperator.EQUAL

            self.constraints[i] = (lhs, comparison, rhs)

    def get_standard_form(self):
        new_LP = copy.deepcopy(self)
        new_LP.add_artificial_variables()
        return new_LP

    def get_TeX(self):
        TeX = []
        TeX.append(f'\\begin{{quote}}')

        min_or_max_str = ''

        if self.min_or_max == MinOrMax.MIN:
            min_or_max_str = 'Minimize'
        elif self.min_or_max == MinOrMax.MAX:
            min_or_max_str = 'Maximize'

        TeX.append(f'  {min_or_max_str} ${TeXify_linear_combination(self.objective)}$ \\\\')
        TeX.append(f'  s.t.')
        TeX.append(f'  \\begin{{align*}}')

        for (lhs, comparison, rhs) in self.constraints:
            if comparison == ComparisonOperator.EQUAL:
                cmp_str = '='
            elif comparison == ComparisonOperator.LESS_EQUAL:
                cmp_str = '\\leq'
            elif comparison == ComparisonOperator.GREATER_EQUAL:
                cmp_str = '\\geq'

            TeX.append(f'    {TeXify_linear_combination(lhs)} &{cmp_str} {rhs.to_TeX()}, \\\\')

        if not self.is_dual:
            TeX.append(f'    {", ".join(self.variable_names)} &\\geq 0.')        
        
        TeX.append(f'  \\end{{align*}}')
        TeX.append(f'\\end{{quote}}')

        return '\n'.join(TeX)

    def get_dual(self):
        if not self.standard_form:
            return self.get_standard_form().get_dual()

        dual = LinearProgram()
        dual.is_dual = True
        dual.min_or_max = MinOrMax.MAX
        dual.variable_names = ['p_' + str(i + 1) for i in range(len(self.constraints))]

        objective = []
        constraints = []
        rows = []
        
        for i, (lhs, comparison, rhs) in enumerate(self.constraints):
            objective.append((Q(rhs), dual.variable_names[i]))

            row = [Q(0) for i in range(len(self.variable_names))]

            for (coef, variable) in lhs:
                row[self.variable_names.index(variable)] = coef

            rows.append(row)

        A_T = Matrix(*rows).transpose()
        c = [Q(0) for i in range(A_T.m)]

        for (coef, variable) in self.objective:
            c[self.variable_names.index(variable)] = coef

        for i in range(A_T.m):
            lhs = []
            comparison = ComparisonOperator.LESS_EQUAL
            rhs = c[i]

            for j in range(A_T.n):
                if A_T.at(i, j).is_nonzero():
                    lhs.append((A_T.at(i, j), dual.variable_names[j]))

            constraints.append((lhs, comparison, rhs))

        dual.objective = objective
        dual.constraints = constraints

        return dual

    def solve(self, printTeX = True):
        LP_standard = self.get_standard_form()
        initial_basis = find_initial_basic_solution(LP_standard, printTeX = printTeX)
        return simplex(LP_standard, initial_basis, SimplexAlg.PRIMAL, printTeX = printTeX)
