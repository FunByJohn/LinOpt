from linear_program import *
from integer_program import *

'''
LP = LinearProgram()
LP.set_objective("min x_2 + x_3")
LP.add_constraint("x_1 + x_2 + x_3 >= 4")
LP.add_constraint("4x_1 + 2x_2 + x_3 <= 12")

solution, basis, tableau = LP.solve(printTeX = False)

print(solution)

# In Mathematica:
# Minimize[{x2 + x3, x1 + x2 + x3 >= 4 && 4 x1 + 2 x2 + x3 <= 12 && x1 >= 0 && x2 >= 0 && x3 >= 0}, {x1, x2, x3}]

IP = IntegerProgram()
IP.set_objective("min x_2 + x_3")
IP.add_constraint("x_1 + x_2 + x_3 >= 4")
IP.add_constraint("4x_1 + 2x_2 + x_3 <= 12")

solution = randomized_gomory(IP, printTeX = False)

print(solution)

# In Mathematica:
# Minimize[{x2 + x3, 
#  x1 + x2 + x3 >= 4 && 4 x1 + 2 x2 + x3 <= 12 && x1 >= 0 && x2 >= 0 &&
#    x3 >= 0 && x1 \[Element] Integers && x2 \[Element] Integers && 
#   x3 \[Element] Integers}, {x1, x2, x3}]
'''

