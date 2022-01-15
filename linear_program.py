from enum import Enum, auto
from rational import *

import re
import copy

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

    def set_objective(self, string):
        (min_max, lin_comb) = string.lstrip().rstrip().split(' ', maxsplit = 1)

        if min_max == 'min':
            self.min_or_max = MinOrMax.MIN
        elif min_max == 'max':
            self.min_or_max = MinOrMax.MAX
        else:
            print('Error when initializing LinearProgram: "min" or "max" not supplied in objective function')
            assert False

        self.objective = self.parse_linear_combination(lin_comb)

    def add_constraint(self, string):
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



        pass

    def add_artificial_variables(self):
        self.standard_form = True

        def get_new_variable_name():
            last_variable_index = int(self.variable_names[-1].split('_')[1])
            new_variable = f'x_{last_variable_index + 1}'
            self.ensure_variable_exists(new_variable)
            return new_variable

        for i in range(len(self.constraints)):
            (lhs, comparison, rhs) = self.constraints[i]

            if comparison == ComparisonOperator.LESS_EQUAL:
                # add slack variable
                lhs.append((Q(1), get_new_variable_name()))
            elif comparison == ComparisonOperator.GREATER_EQUAL:
                # add surplus variable
                lhs.append((Q(-1), get_new_variable_name()))

            comparison = ComparisonOperator.EQUAL

            self.constraints[i] = (lhs, comparison, rhs)

    def get_standard_form(self):
        new_LP = copy.deepcopy(self)
        new_LP.add_artificial_variables()
        return new_LP

    def get_TeX(self):
        
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

        TeX.append(f'    {", ".join(self.variable_names)} &\\geq 0.')        
        TeX.append(f'  \\end{{align*}}')
        TeX.append(f'\\end{{quote}}')

        return '\n'.join(TeX)

    def get_dual(self):
        pass


LP = LinearProgram()
LP.set_objective('  min  2 x_1   + x_2 + 2 x_3')
LP.add_constraint('        x_1         + 3 x_3  <= 5 ')
LP.add_constraint('              2 x_2 +   x_3  <= 3 ')
LP.add_constraint('      2 x_1 +   x_2 +   x_3  >= 2 ')
#LP.add_constraint('        x_1                  >= 0 ')
#LP.add_constraint('                x_2          >= 0 ') # TODO: How to handle x_1, x_2, x_3 >= 0 constraint?
#LP.add_constraint('                        x_3  >= 0 ')

LP_standard = LP.get_standard_form()
print(LP_standard.get_TeX())