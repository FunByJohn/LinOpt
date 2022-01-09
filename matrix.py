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
        assert scalar.is_nonzero()

        scale = Matrix_identity(self.m)
        scale.entries[i][i] = Q(scalar)

        return Matrix_mul(scale, self)

    def add_rows(self, i, j, scalar):
        add = Matrix_identity(self.m)
        add.entries[i][j] = Q(scalar)

        return Matrix_mul(add, self)

    def rref(self):
        return Matrix_rref(Matrix(*self.entries))

    def inverse(self):
        assert self.m == self.n

        # TODO: Check if the matrix was invertible

        return Matrix_rref(Matrix_combine_horizontal(self, Matrix_identity(self.n))) \
               .submatrix(self.n, 0, self.n, self.n)

    def submatrix(self, left, top, width, height):
        rows = []

        for i in range(top, top + height):
            row = []
            for j in range(left, left + width):
                row.append(self.at(i, j))
            rows.append(row)

        return Matrix(*rows)


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

def Matrix_rref(matrix):
    left = 0
    top = 0
    pivots = []
    
    while left < matrix.n and top < matrix.m:
        nonzero_index = None
        for i in range(top, matrix.m):
            if matrix.at(i, left).is_nonzero():
                nonzero_index = i
                break
        
        if nonzero_index == None:
            left += 1
        else:
            pivots.append((top, left))
            matrix = matrix.swap_rows(top, nonzero_index)
            matrix = matrix.scale_row(top, matrix.at(top, left).mul_inv())

            for i in range(top + 1, matrix.m):
                matrix = matrix.add_rows(i, top, matrix.at(i, left).add_inv())

            left += 1
            top += 1

    for (pivot_i, pivot_j) in pivots:
        for i in range(0, pivot_i):
            matrix = matrix.add_rows(i, pivot_i, matrix.at(i, pivot_j).add_inv())
    
    return matrix

def Matrix_combine_horizontal(a, b):
    assert a.m == b.m

    rows = []
    for i in range(a.m):
        row = []
        for j in range(a.n + b.n):
            if j < a.n:
                row.append(a.at(i, j))
            else:
                row.append(b.at(i, j - a.n))
        rows.append(row)

    return Matrix(*rows)

def Matrix_combine_vertical(a, b):
    # I feel like a criminal for doing this
    return Matrix_combine_horizontal(a.transpose(), b.transpose()).transpose()

