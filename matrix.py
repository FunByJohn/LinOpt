from rational import *

class Matrix:

    def __init__(self, *rows):
        self.entries = []

        for line in rows:
            row = []

            for entry in line:
                row.append(Q(entry))
            
            self.entries.append(row)

        self.m = len(self.entries)
        self.n = len(self.entries[0])

    def __str__(self):
        return str(self.entries)

    def at(self, i, j):
        return self.entries[i][j]

    def transpose(self):
        rows = []

        for i in range(self.n):
            row = []
            for j in range(self.m):
                row.append(self.at(j, i))
            rows.append(row)

        return Matrix(*rows)

    def swap_rows(self, i, j):
        swap = Matrix_identity(self.m)
        swap.entries[i][i] = Q(0)
        swap.entries[j][j] = Q(0)
        swap.entries[i][j] = Q(1)
        swap.entries[j][i] = Q(1)

        return Matrix_mul(swap, self)

    def scale_row(self, i, scalar):
        assert scalar.p != 0

        scale = Matrix_identity(self.m)
        scale.entries[i][i] = Q(scalar)

        return Matrix_mul(scale, self)

    def add_rows(self, i, j, scalar):
        add = Matrix_identity(self.m)
        add.entries[i][j] = Q(scalar)

        return Matrix_mul(add, self)

    #def rref(self):
    #    def alg(matrix):
    #        left = 0
    #        top = 0
    #
    #        while left != matrix.n:
    #            for i in range(top, matrix.m):
    #                #
    #
    #        return matrix
    #
    #    return alg(Matrix(*self.entries))

def Matrix_identity(n):
    rows = []

    for i in range(n):
        row = []
        for j in range(n):
            row.append(1 if i == j else 0)
        rows.append(row)

    return Matrix(*rows)

def Matrix_mul(a, b):
    assert a.n == b.m

    rows = []
    for i in range(a.m):
        row = []
        for j in range(b.n):
            accum = Q(0)
            for k in range(a.n):
                accum = Q_add(accum, Q_mul(a.at(i, k), b.at(k, j)))
            row.append(accum)
        rows.append(row)

    return Matrix(*rows)