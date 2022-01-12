import math

from collections import deque
from graph import *

def bellman_ford(G, vertex_name, printTeX = True):
    n = G.vert_names.index(vertex_name)

    # Set weights that are 0 to infinity to make the algorithm in the book work
    for i in range(len(G.vert_names)):
        for j in range(len(G.vert_names)):
            if (G.weight_matrix[i][j] == 0):
                G.weight_matrix[i][j] = math.inf

    # p[i][j] is the book's p_i(j)
    p = [[math.inf for i in range(len(G.vert_names))] for j in range(len(G.vert_names))]

    # This is s^* from Kent's note on Bellman-Ford
    s = [None for i in range(len(G.vert_names))]

    # Set p_n(t) = 0 for all t
    for t in range(len(G.vert_names)):
        p[n][t] = 0

    # Iterate
    for t in range(1, len(G.vert_names)):
        for i in range(len(G.vert_names)):
            if i != n:
                value = math.inf
                through = None

                for k in range(len(G.vert_names)):
                    if G.weight_matrix[i][k] != math.inf:
                        new_value = G.weight_matrix[i][k] + p[k][t - 1]

                        if new_value <= value:
                            value = new_value
                            through = k

                p[i][t] = value

                if p[i][t] != p[i][t-1]:
                    s[i] = through

    # Compute paths
    paths = [[] for i in range(len(G.vert_names))]

    for i in range(len(G.vert_names)):
        k = i
        while (s[k] != None):
            paths[i].append(G.vert_names[k])
            k = s[k]

        paths[i].append(G.vert_names[n])

    return paths

def dijkstra(G, vertex_name, printTeX = True):
    '''
    Implementation of Dijkstras algorithm as described on p. 341, but using bellman_ford to actually print out the paths
    '''
    TeX = []
    active_indices = list(map(lambda x: G.vert_names.index(x), G.vert_names))
    n = G.vert_names.index(vertex_name)
    path_lengths = [0 for i in range(len(active_indices))]
    active_indices.remove(n)
    parent = [None for i in range(len(G.vert_names))]

    def TeX_repr(n):
        if n == math.inf:
            return '\\infty'
        return n

    if printTeX:
        TeX.append(f'\\begin{{itemize}}')
        TeX.append(f'  \\item Vi følger algoritmen på s. 341 i bogen. Lad $n = {G.vert_names[n]}$.')

    # Set weights that are 0 to infinity to make the algorithm in the book work
    for i in range(len(G.vert_names)):
        for j in range(len(G.vert_names)):
            if (G.weight_matrix[i][j] == 0):
                G.weight_matrix[i][j] = math.inf

    # Bellman-Ford to find the actual paths, since the book doesn't specify how to find these in a straight forward way
    paths = bellman_ford(G, vertex_name, False)

    while (len(active_indices) > 0):
        if printTeX:
            TeX.append(f'  \\item De knuder vi mangler at betragte er $\\{{{", ".join(list(map(lambda x: G.vert_names[x], active_indices)))}\\}}$.')
        
        # Step 1: Find a node l ≠ n such that weights[l][n] <= weights[i][n] for all i != n. Set p^*_l = weights[l][n].
        l = -1
        l_weight = math.inf

        for i in active_indices:
            if G.weight_matrix[i][n] <= l_weight:
                l = i
                l_weight = G.weight_matrix[i][n]

        if printTeX:
            TeX.append(f'  \\item \\textbf{{Skridt 1.}} Vi har at $\\ell = {G.vert_names[l]}$, og får at $p_{{{G.vert_names[l]}}}^* = {G.weight_matrix[l][n]}$.')

        path_lengths[l] = G.weight_matrix[l][n]

        # Step 2: For every node i != l, n, set weights[i][n] = min(weights[i][n], weights[i][l] + weights[l][n])
        step2TeX = []

        for i in active_indices:
            if i != l:
                #if G.weight_matrix[i][l] + G.weight_matrix[l][n] <= G.weight_matrix[i][n]:
                #    paths[i] = [l] + paths[i]

                if printTeX:
                    step2TeX.append(f'c_{{{G.vert_names[i]}{G.vert_names[n]}}} &:= \\text{{min}}' + \
                                    f'(c_{{{G.vert_names[i]}{G.vert_names[n]}}}, c_{{{G.vert_names[i]}{G.vert_names[l]}}} + ' + \
                                    f'c_{{{G.vert_names[l]}{G.vert_names[n]}}}) = \\text{{min}}({TeX_repr(G.weight_matrix[i][n])}, ' + \
                                    f'{TeX_repr(G.weight_matrix[i][l] + G.weight_matrix[l][n])}) = ' + \
                                    f'{TeX_repr(min(G.weight_matrix[i][n], G.weight_matrix[i][l] + G.weight_matrix[l][n]))}.')

                G.weight_matrix[i][n] = min(G.weight_matrix[i][n], G.weight_matrix[i][l] + G.weight_matrix[l][n])

        if printTeX and len(step2TeX) > 0:
            split_str = ' \\\\ \n'
            TeX.append(f'  \\item \\textbf{{Skridt 2.}} Vi sætter \\begin{{align*}} {split_str.join(step2TeX)} \\end{{align*}}')

        # Step 3: Remove node l from the graph and apply the same steps to the new graph.
        for i in active_indices:
            G.weight_matrix[i][l] = math.inf
            G.weight_matrix[l][i] = math.inf

        if printTeX:
            TeX.append(f'  \\item \\textbf{{Skridt 3.}} Vi fjerner nu knuden ${G.vert_names[l]}$ ved at sætte $c_{{i{G.vert_names[l]}}} = c_{{{G.vert_names[l]}i}} = \\infty$ for alle $i$.')

        active_indices.remove(l)

    if printTeX:
        TeX.append(f'  \\item Der er nu ikke flere knuder at betragte, hvormed algoritmen terminerer.')

        for i in range(len(G.vert_names)):
            path = paths[i]
            split_str = " \\to "

            if i != n:
                label = G.vert_names[i]
                TeX.append(f'  \\item Korteste vej fra ${G.vert_names[i]}$ til ${G.vert_names[n]}$ er ${split_str.join(path)}$ med længden $p_{{{label}}}^* = {path_lengths[i]}$.')

        TeX.append(f'\\end{{itemize}}')

        print('\n'.join(TeX))

G = Graph(
    '1 2 3 4 5 6',
    '(1, 2)=1, (1, 3)=1, (1, 4)=3, (2, 5)=4, (3, 2)=2, (3, 4)=4, (3, 5)=5, (4, 6)=9, (5, 4)=1, (5, 6)=6',
)

dijkstra(G, '6', True)