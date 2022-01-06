from math import gcd

class Q:

    def __init__(self, p, q = 1):
        assert q != 0

        g = gcd(p, q)

        self.p = p // g
        self.q = q // g

    def __str__(self):
        return self.tostring()

    def __repr__(self):
        return self.tostring()

    def add_inv(self):
        return Q(-self.p, self.q)

    def mul_inv(self):
        return Q(self.q, self.p)

    def tostring(self):
        if self.q == 1:
            return str(self.p)

        return f'{self.p}/{self.q}'

def Q_add(lhs, rhs):
    return Q(rhs.p * lhs.q + lhs.p * rhs.q, lhs.q * rhs.q)

def Q_sub(lhs, rhs):
    return Q_add(lhs, rhs.add_inv())

def Q_mul(lhs, rhs):
    return Q(lhs.p * rhs.p, lhs.q * rhs.q)

def Q_div(lhs, rhs):
    return Q_mul(lhs, rhs.mul_inv())