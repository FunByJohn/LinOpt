import re

class Graph:

    def __init__(self, verts_string, edges_string):
        self.directed = True
        self.edge_pairs = []
        self.edge_weights = []

        self.setup_verts(verts_string)
        self.setup_edges(edges_string)
        self.setup_matrices()

    def setup_verts(self, verts_string):
        self.vert_names = re.split(r'[\s,]+', verts_string.lstrip().rstrip())

    def setup_edges(self, edges_string):
        edges = re.findall(r'\(\s*(\w+)\s*,\s*(\w+)\s*\)\s*=\s*([-]{0,1}\s*\d+)', edges_string)

        for (src, dst, weight) in edges:
            src = self.vert_names.index(src)
            dst = self.vert_names.index(dst)
            weight = int(weight)
            
            self.edge_pairs.append((src, dst))
            self.edge_weights.append((src, dst, weight))

    def setup_matrices(self):
        def adj_matrix_entry(i, j):
            return 1 if (i, j) in self.edge_pairs else 0

        self.adj_matrix = [[adj_matrix_entry(i, j) for j in range(len(self.vert_names))] for i in range(len(self.vert_names))]
        self.weight_matrix = [[0 for j in range(len(self.vert_names))] for i in range(len(self.vert_names))]

        for (src, dst, weight) in self.edge_weights:
            self.weight_matrix[src][dst] = weight

    def make_undirected(self):
        self.directed = False

        edge_pairs_to_add = []
        edge_weights_to_add = []

        for (src, dst) in self.edge_pairs:
            edge_pairs_to_add.append((dst, src))

        for (src, dst, weight) in self.edge_weights:
            edge_weights_to_add.append((dst, src, weight))

        self.edge_pairs = self.edge_pairs + edge_pairs_to_add
        self.edge_weights = self.edge_weights + edge_weights_to_add
        self.setup_matrices()

    def present(self):
        print(self.vert_names)
        print(self.edge_pairs)
        print(self.adj_matrix)
        print(self.weight_matrix)

