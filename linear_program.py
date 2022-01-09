import re

class LinearProgram:

    def __init__(self):
        self.variable_names = None
        self.constraints = []
        self.slack_surplus_counter = 1

    def set_variables(self, string):
        self.variable_names = re.split(r'\s+', string.lstrip().rstrip())

    def set_objective(self, string):
        pass

    def add_constraint(self, string):
        comparison = re.search(r'([<>]{0,1}=)', string).group(0)
        (lhs, rhs) = re.split(r'\s*[<>]{0,1}=\s*', string)
        lhs = self.parse_linear_combination(lhs)
        rhs = int(rhs)

        if comparison == '<=':
            print('add slack variable')
        else:
            print('add surplus variable')

        print([lhs, comparison, rhs])
        pass

    def parse_linear_combination(self, string):
        result = []
        parts = re.findall(r'[-\+]{0,1}\s*\d*\s*[a-z]_\d+', string)

        for part in parts:
            match = re.search(r'([-\+]{0,1}\s*\d*)\s*([a-z]_\d+)', part)
            coef = match.group(1).lstrip().rstrip()
            name = match.group(2)
            sgn = 1 if re.match(r'-', coef) == None else -1

            if coef == '-':
                coef = -1
            else:
                coef = re.sub(r'[-\+]', '', coef)

                if coef == '':
                    coef = 1
                else:
                    coef = sgn * int(coef)

            result.append((coef, name))

        return result

    def present(self):
        pass

    def dual(self):
        pass


LP = LinearProgram()
LP.set_variables(' x_1 x_2 x_3 ')
LP.set_objective('  min  2 x_1   - x_2 + 2 x_3')
LP.add_constraint('      - x_1         + 3 x_3  <= 5 ')
LP.add_constraint('              2 x_2 +   x_3  <= 5 ')
LP.add_constraint('      2 x_1 +   x_2 +   x_3  >= 2 ')
LP.add_constraint('              4 x_2 - 3 x_3   = 7 ')
LP.add_constraint('        x_1                  >= 0 ')
LP.add_constraint('                x_2          >= 0 ')
LP.add_constraint('                        x_3  >= 0 ')

#LP_dual = LP.dual()
#print(LP_dual.present())