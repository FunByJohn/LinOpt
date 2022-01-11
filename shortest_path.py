import math

from collections import deque
from graph import *

def dijkstra(G, vertex_name):
    '''
    Implementation of Dijkstras algorithm as described on p. 341, with a modification to print out shortest paths
    '''
    active_indices = list(map(lambda x: G.vert_names.index(x), G.vert_names))
    n = G.vert_names.index(vertex_name)
    paths = [[n] for i in range(len(active_indices))]
    path_lengths = [0 for i in range(len(active_indices))]
    active_indices.remove(n)
    discovered_order = deque()

    # Set weights that are 0 to infinity to make the algorithm in the book work
    for i in range(len(G.vert_names)):
        for j in range(len(G.vert_names)):
            if (G.weight_matrix[i][j] == 0):
                G.weight_matrix[i][j] = math.inf

    while (len(active_indices) > 0):
        # Step 1: Find a node l ≠ n such that weights[l][n] <= weights[i][n] for all i != n. Set p^*_l = weights[l][n].
        l = -1
        l_weight = math.inf

        for i in active_indices:
            if G.weight_matrix[i][n] <= l_weight:
                l = i
                l_weight = G.weight_matrix[i][n]

        path_lengths[l] = G.weight_matrix[l][n]

        # Step 2: For every node i != l, n, set weights[i][n] = min(weights[i][n], weights[i][l] + weights[l][n])
        for i in active_indices:
            if i != l:
                G.weight_matrix[i][n] = min(G.weight_matrix[i][n], G.weight_matrix[i][l] + G.weight_matrix[l][n])

        # Step 3: Remove node l from the graph and apply the same steps to the new graph.
        for i in active_indices:
            G.weight_matrix[i][l] = math.inf
            G.weight_matrix[l][i] = math.inf

        active_indices.remove(l)

    # TODO: Actually discover the paths

    for i in range(len(G.vert_names)):
        path = paths[i]
        path = list(map(lambda x: G.vert_names[x], path))
        print(f'Shortest path from {G.vert_names[i]} to {G.vert_names[n]}: {path}')

    print(path_lengths)

G = Graph(
    '1 2 3 4 5 6',
    '(1, 2)=1, (1, 3)=1, (1, 4)=3, (2, 5)=4, (3, 2)=2, (3, 4)=4, (3, 5)=5, (4, 6)=9, (5, 4)=1, (5, 6)=6',
)

dijkstra(G, '6')