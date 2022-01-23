from matrix import *
from graph import *

class Network(Graph):

    def __init__(self, verts_string, edges_string):
        self.directed = True
        self.edge_pairs = []
        self.edge_weights = []
        self.vertex_weights = []

        if '(' in verts_string:
            self.vert_names = []
            
            for (vert, weight) in re.findall(r'\(\s*(\w+)\s*\)\s*=\s*([-]{0,1}\s*\d+)', verts_string):
                self.vert_names.append(vert)
                self.vertex_weights.append(int(weight))
        else:
            self.setup_verts(verts_string)

        self.setup_edges(edges_string)
        self.setup_matrices()

    def make_tree(self, edges_string):
        tree = []

        for (lhs, rhs) in re.findall(r'\(\s*(\w+)\s*,\s*(\w+)\s*\)', edges_string):
            tree.append((self.vert_names.index(lhs), self.vert_names.index(rhs)))

        return tree

    def make_A(self):
        A = Matrix(*[[Q(0) for j in range(len(self.edge_pairs))] for i in range(len(self.vert_names))])

        for i in range(len(self.vert_names)):
            for k in range(len(self.edge_pairs)):
                (start_node, end_node) = self.edge_pairs[k]

                if i == start_node:
                    A.entries[i][k] = Q(1)
                elif i == end_node:
                    A.entries[i][k] = Q(-1)
        
        return A

def network_simplex(network, tree, printTeX = True):
    # We follow the algorithm on p. 286

    if printTeX:
        tree_str = []

        for (node_start, node_end) in tree:
            tree_str.append(f'({network.vert_names[node_start]}, {network.vert_names[node_end]})')

        tree_str = ', '.join(tree_str)
        print(f'We now have the tree ${tree_str}$ \\\\')

    # Step 1
    A = network.make_A()
    A_trunc = A.submatrix(0, 0, A.n, A.m - 1)
    A_trunc_square = None
    b_trunc = Matrix(*[[Q(network.vertex_weights[i])] for i in range(len(network.vertex_weights) - 1)])
    tree_indices = []

    for edge in tree:
        index = network.edge_pairs.index(edge)
        tree_indices.append(index)
        column = A_trunc.submatrix(index, 0, 1, A_trunc.m)

        if A_trunc_square == None:
            A_trunc_square = column
        else:
            A_trunc_square = Matrix_combine_horizontal(A_trunc_square, column)

    solved_flow_values = Matrix_mul(A_trunc_square.inverse(), b_trunc)
    f = [0 for i in range(len(network.edge_pairs))]

    for i in range(len(tree_indices)):
        f[tree_indices[i]] = solved_flow_values.at(i, 0).to_int()

    def set_f(i, j, value):
        f[network.edge_pairs.index((i, j))] = value

    def get_f(i, j):
        return f[network.edge_pairs.index((i, j))]

    if printTeX:
        string = []
        for (node_start, node_end) in network.edge_pairs:
            string.append(f'$f_{{{network.vert_names[node_start]},{network.vert_names[node_end]}}} = {get_f(node_start, node_end)}$')
        print('This gives us the basic feasible solution \\\\')
        print(', '.join(string) + ' \\\\')

    # Step 2
    if printTeX:
        print('Now we compute the dual vector. The equation system is \\\\')
    
    n = len(network.vert_names)
    M = Matrix(*[[Q(0) for i in range(n)] for j in range(n)])
    v = Matrix(*[[Q(0)] for i in range(n)])

    for row, (i, j) in enumerate(tree):
        M.entries[row][i] = Q(1)
        M.entries[row][j] = Q(-1)
        v.entries[row][0] = Q(network.weight_matrix[i][j])

        if printTeX:
            print(f'$p_{network.vert_names[i]} - p_{network.vert_names[j]} = {network.weight_matrix[i][j]}$, \\\\')

    M.entries[n - 1][n - 1] = Q(1)

    if printTeX:
        print(f'$p_{network.vert_names[n - 1]} = 0$. \\\\')

    p = Matrix_mul(M.inverse(), v)
    p = [p.at(i, 0).to_int() for i in range(n)]
    
    if printTeX:
        p_str = ', '.join([str(x) for x in p])
        print(f'This has the solution $p = ({p_str})$. \\\\')
    
    # Step 3
    reduced_costs = [[None for i in range(n)] for j in range(n)]
    reduced_costs_to_check = []
    
    for i in range(n):
        for j in range(n):
            if (i, j) in network.edge_pairs:
                reduced_costs[i][j] = network.weight_matrix[i][j] - (p[i] - p[j])

                if (i, j) in tree:
                    reduced_costs_to_check.append(((i, j, reduced_costs[i][j])))

    if printTeX:
        print("The reduced costs $\\overline{{c}}_{{ij}}$ are: \\\\")

        lines = []
        for i in range(n):
            line = []
            for j in range(n):
                if reduced_costs[i][j] != None:
                    line.append(f'$\\overline{{c}}_{{{network.vert_names[i]},{network.vert_names[j]}}} = {reduced_costs[i][j]}$')
            
            if len(line) > 0:
                lines.append(', '.join(line))

        print(' \\\\ \n'.join(lines) + ' \\\\')

    is_optimal = True
    for (i, j, cost) in reduced_costs_to_check:
        pass

    if is_optimal:
        if printTeX:
            print('All reduced costs outside our tree are nonnegative, which means that the basic feasible solution $f_{{ij}}$ is optimal.')

    return ":O"

def ford_fulkerson(network):
    pass

# Network simplex
N = Network(
    '(1)=8, (2)=-2, (3)=5, (4)=-7, (5)=-4',
    '(1,2)=2, (1,4)=9, (2,3)=1, (2,5)=4, (3,1)=6, (3,4)=8, (3,5)=5, (5,1)=-7, (5,4)=3'
)

T = N.make_tree('(1, 2), (3, 1), (1, 4), (3, 5)')

print(network_simplex(N, T))

# Ford Fulkerson
N = Network(
    'A B C D E F G',
    '(A,B)=5, (A,C)=7, (A,D)=4, (B,C)=1, (B,E)=3, (C,D)=2, (C,E)=4, (C,F)=5, (D,F)=4, (E,F)=1, (E,G)=9, (F,G)=6'
)

print(ford_fulkerson(N))