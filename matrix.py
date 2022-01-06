import sys
import re
from rational import *

class Matrix:

    def __init__(self, *rows):
        self.entries = []

        for line in rows:
            row = []

            for entry in line:
                if isinstance(entry, Q):
                    row.append(Q(entry.p, entry.q))
                elif isinstance(entry, int):
                    row.append(Q(entry, 1))
                elif isinstance(entry, str):
                    if '/' in entry:
                        (p, q) = re.split('\s*/\s*', entry)
                        row.append(Q(int(p), int(q)))
                    else:
                        row.append(Q(int(entry)))
                else:
                    raise SystemExit('Error in matrix initialization: unsupported type given as input')
            
            self.entries.append(row)

        self.m = len(self.entries)
        self.n = len(self.entries[0])

    def at(self, i, j):
        return self.entries[i][j]

    def __str__(self):
        return str(self.entries)

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