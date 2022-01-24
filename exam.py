from linear_program import *
from integer_program import *
from simplex import *
from flow import *

############################################################
# Opg 1
############################################################

# (a)

LP = LinearProgram()
LP.set_objective("max 2x_1 + 3x_2 + 4x_3")
LP.add_constraint("x_1 + x_2 + x_3 = 5")
LP.add_constraint("2x_1 + x_2 + 2x_3 >= 5")
LP.add_constraint("x_1 + 2x_2 + 2x_3 <= 8") # Actually (8-t)

print(LP.get_dual().get_TeX())

# (b)

initial_basis = find_initial_basic_solution(LP, printTeX = True)

# (c)

solution, current_basis, tableau = simplex(LP.get_standard_form(), initial_basis, SimplexAlg.PRIMAL, printTeX = True)

# (g)

LP = LinearProgram()
LP.set_objective("max 2x_1 + 3x_2 + 4x_3")
LP.add_constraint("x_1 + x_2 + x_3 = 5")
LP.add_constraint("2x_1 + x_2 + 2x_3 >= 5")
LP.add_constraint("x_1 + 2x_2 + 2x_3 <= 8")
LP_standard = LP.get_standard_form()

print(LP_standard.get_TeX())

LP = LinearProgram()
LP.set_objective("max 2x_1 + x_2 + x_3 + 3x_6")
LP.add_constraint("x_1 + x_2 + x_3 + x_6 = 5")
LP.add_constraint("2x_1 + x_2 + 2x_3 - x_4 + x_6 = 5")
LP.add_constraint("x_1 + 2x_2 + 2x_3 + x_5 + x_6 = 8")

simplex(LP, current_basis, SimplexAlg.PRIMAL, printTeX = True)

############################################################
# Opg 2
############################################################

LP = LinearProgram()
LP.set_objective("max x_2 + x_3 + x_4")
LP.add_constraint("x_1 + x_2 + x_5 <= 1")
LP.add_constraint("x_2 + x_3 = 1")
LP.add_constraint("x_3 + x_4 = 1")
LP.add_constraint("x_2 + x_4 + x_5 = 1")

# (a)

print("Is feasible?")
print(LP.is_feasible([Q("0"), Q("1/2"), Q("1/2"), Q("1/2"), Q("0")]))

LP_standard = LP.get_standard_form()

basis = [1, 2, 3, 5]

simplex(LP_standard, basis, SimplexAlg.PRIMAL, printTeX = True)

# (b)

HP = IntegerProgram()
HP.set_objective("max x_2 + x_3 + x_4")
HP.add_constraint("x_1 + x_2 + x_5 <= 1")
HP.add_constraint("x_2 + x_3 = 1")
HP.add_constraint("x_3 + x_4 = 1")
HP.add_constraint("x_2 + x_4 + x_5 = 1")

while True:
    try:
        print(randomized_gomory(HP, printTeX = True))
        break
    except:
        print("===")
        print("Runtime error in randomized_gomory, we try again!")
        print("===")
        continue

# (d)

LP_standard.add_constraint("x_2 = 0")
#LP_standard.add_constraint("x_3 - x_5 = 0")
#LP_standard.add_constraint("x_4 = 0")
#LP_standard.add_constraint("x_1 + x_6 = 0")

print(LP_standard.get_TeX())
print(LP_standard.solve())

############################################################
# Opg 3
############################################################

N = Network(
    '(s1)=36, (s2)=42, (s3)=42, (d1)=-45, (d2)=-23, (d3)=-36, (d4)=-16',
    '(s1, d1) = 6, (s1, d2) = 5, (s1, d3) = 11, (s1, d4) = 4,' + \
    '(s2, d1) = 2, (s2, d2) = 8, (s2, d3) = 9, (s2, d4) = 6,' + \
    '(s3, d1) = 10, (s3, d2) = 3, (s3, d3) = 8, (s3, d4) = 2'
)

T = N.make_tree('(s1,d1), (s1,d4), (s2,d1), (s2,d3), (s3,d2), (s3,d3)')

cycles = [
    [(0,4), (2,4), (2,5), (1,5), (1,3), (0,3)]
]

Fs = [
    [(0,4),(2,5),(1,3)]
]

Bs = [
    [(2,4),(1,5),(0,3)]
]

print(network_simplex(N, T, cycles, Fs, Bs, printTeX = True))

N = Network(
    '(s1)=36, (s2)=42, (s3)=42, (d1)=-45, (d2)=-23, (d3)=-36, (d4)=-16',
    '(s1, d1) = 6, (s1, d2) = 5, (s1, d3) = 11, (s1, d4) = 4,' + \
    '(s2, d1) = 2, (s2, d2) = 8, (s2, d3) = 9, (s2, d4) = 6,' + \
    '(s3, d1) = 10, (s3, d2) = 3, (s3, d3) = 8, (s3, d4) = 2'
)

T = N.make_tree('(s1,d1), (s2,d1), (s2,d2), (s2,d3), (s3,d3), (s3,d4)')

cycles = [
    
]

Fs = [
    
]

Bs = [
    
]

print(network_simplex(N, T, cycles, Fs, Bs, printTeX = True))