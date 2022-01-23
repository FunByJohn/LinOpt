import re

from math import gcd

class Q:

    def __init__(self, p, q = 1):
        if isinstance(p, Q):
            inst = p
            p = inst.p
            q = inst.q
        elif isinstance(p, str):
            if '/' in p:
                (p_str, q_str) = re.split('\s*/\s*', p)
                p = int(p_str)
                q = int(q_str)
            else:
                p = int(p)
        elif not isinstance(p, int):
            print('Error when initializing rational: unsupported type given as input')
            assert False

        assert q != 0

        g = gcd(p, q)

        self.p = p // g
        self.q = q // g

        if self.q < 0:
            self.p = -self.p
            self.q = -self.q

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return self.to_string()

    def add_inv(self):
        return Q(-self.p, self.q)

    def mul_inv(self):
        return Q(self.q, self.p)

    def to_string(self):
        if self.p == 0:
            return 'Q(0)'

        if self.q == 1:
            return f'Q({str(self.p)})'

        return f'Q({self.p}/{self.q})'

    def to_TeX(self):
        if self.p == 0:
            return '0'

        if self.q == 1:
            return str(self.p)

        return f'{"-" if self.is_negative() else ""}\\tfrac{{{abs(self.p)}}}{{{self.q}}}'

    def to_int(self):
        assert self.q == 1
        return self.p

    def is_zero(self):
        return self.p == 0

    def is_one(self):
        return self.p == 1 and self.q == 1

    def is_nonzero(self):
        return self.p != 0

    def is_negative(self):
        return self.p < 0

    def is_positive(self):
        return self.p > 0

    def is_integer(self):
        return self.q == 1

def Q_add(lhs, rhs):
    return Q(rhs.p * lhs.q + lhs.p * rhs.q, lhs.q * rhs.q)

def Q_sub(lhs, rhs):
    return Q_add(lhs, rhs.add_inv())

def Q_mul(lhs, rhs):
    return Q(lhs.p * rhs.p, lhs.q * rhs.q)

def Q_div(lhs, rhs):
    return Q_mul(lhs, rhs.mul_inv())

def Q_equals(lhs, rhs):
    return lhs.p == rhs.p and lhs.q == rhs.q

def Q_le(lhs, rhs):
    return lhs.p * rhs.q < rhs.p * lhs.q

def Q_ge(lhs, rhs):
    return lhs.p * rhs.q > rhs.p * lhs.q

def Q_leq(lhs, rhs):
    return Q_equals(lhs, rhs) or Q_le(lhs, rhs)

def Q_geq(lhs, rhs):
    return Q_equals(lhs, rhs) or Q_ge(lhs, rhs)

def Q_abs(x):
    return Q(-x.p, x.q) if x.p < 0 else Q(x)

def Q_floor(x):
    if x.is_negative():
        return Q_sub(Q_floor(x.add_inv()).add_inv(), Q(1))
    
    return Q(x.p // x.q)
    