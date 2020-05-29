import numpy as np
import random
import warnings
from collections.abc import Iterable

DTYPE = np.int64

class PolyRingModP:

    def __init__(self, N, p=2048, coefs=None):
        self.p = p
        # the degree of poly
        self.N = N
        # the coefficients the the polynomial
        # self.coefs[i] indicate the coefficient of x^i
        if type(coefs) is np.ndarray:
            coefs = PolyRingModP.remove_right_zeros(coefs)
            self.coefs = np.pad(coefs, (0, N-coefs.size))
        elif isinstance(coefs, list):
            coefs = PolyRingModP.remove_right_zeros(coefs)
            self.coefs = np.pad(np.array(coefs, dtype=DTYPE), (0, N-len(coefs)))
        elif not coefs:
            self.coefs = np.zeros(N, dtype=DTYPE)

    def __repr__(self):
        return f"{type(self).__name__}({str(self.N)}, {str(self.p)}, {str(self.coefs)})"

    def mod_p(self):
        self.coefs = self.coefs % self.p

    def change_mod(self, q):
        return PolyRingModP(self.N, q, self.coefs)

    # get degree of poly a
    @staticmethod
    def get_degree(a):
        if isinstance(a, Iterable):
            n = len(a)
            for i in range(n-1, -1, -1):
                if a[i]:
                    return i
        else:
            raise NotImplementedError

        return 0

    @staticmethod
    def remove_right_zeros(a):
        if isinstance(a, Iterable):
            return a[:PolyRingModP.get_degree(a)+1]
        else:
            raise NotImplementedError

    @staticmethod
    # extended euclidean algorithm
    def exgcd_not_rec(a, b):
        stack = []
        while b != 0:
            stack.append((a, b))
            a, b = b, a % b
        x = 1
        y = 0
        while (len(stack)):
            a, b = stack.pop()
            x, y = y, x - a // b * y
        return x, y

    # this should be used when d is large
    def random_poly_shuffle(self, d):
        for i in range(d):
            self.coefs[i] = 1
        for i in range(d, 2*d):
            self.coefs[i] = -1
        np.random.shuffle(self.coefs)

    # this should be used if d is small
    def random_poly(self, d):
        """
        randomly choose 1 and -1
        Args:
            d: the number of 1 or -1
        """

        warnings.warn("This is not a standard method!")

        index_num = 2 * d
        assert index_num < self.N
        index = []
        while len(index) <= index_num:
            random_num = random.randint(0, self.N-1)
            if random_num not in index:
                index.append(random_num)

        self.index_p1 = index[:d]
        self.index_n1 = index[d:]

        for i in range(d):
            self.coefs[self.index_p1[i]] = 1
        for i in range(d):
            self.coefs[self.index_n1[i]] = -1

    # if poly only contains {1, -1} and we know the index
    # the multiplication can be simplify
    def fast_mul(self, other):
        """
        This can be only used if coefficients of the two polynomials are 1 or -1
        :param other:
        :return:
        """
        if isinstance(other, int):
            c = np.zeros(self.N, dtype=DTYPE)
            for p_1, n_1 in zip(self.index_p1, self.index_n1):
                c[p_1] = 1
                c[n_1] = -1
            return PolyRingModP(self.N, self.p, coefs=c*other)
        elif isinstance(other, PolyRingModP):
            assert self.N == other.N
            assert self.p == other.p

            c = np.zeros(self.N, dtype=DTYPE)
            index1 = self.index_p1 + self.index_n1
            index2 = other.index_p1 + other.index_n1
            for coef1 in index1:
                for coef2 in index2:
                    if coef1+coef2 >= self.N:
                        c[coef1+coef2-self.N] += self.coefs[coef1] * other.coefs[coef2]
                    else:
                        c[coef1+coef2] += self.coefs[coef1] * other.coefs[coef2]
            return PolyRingModP(self.N, self.p, coefs=c)
        else:
            raise NotImplementedError

    @staticmethod
    def mul(N, a, b):
        """
            calculate the multiplication of two poly
        :param N: X^N-1
        :param p: modular
        :param a: poly1
        :param b: poly2
        :return: a*b \in Z_p[X]/(X^N – 1)
        """
        assert isinstance(a, np.ndarray)
        assert isinstance(b, np.ndarray)
        c = np.zeros(N, dtype=DTYPE)
        for i in range(len(a)):
            for j in range(len(b)):
                if i + j >= N:
                    # maybe overflow
                    c[i + j - N] += a[i] * b[j]
                else:
                    c[i + j] += a[i] * b[j]
        return c

    @staticmethod
    def product_mul(f1, f2, f3):
        """
             calculate f1*f2+f3
        :param f1: the Polynomial only contains {-1,1}
        :param f2: same as f1
        :param f3: same as f1
        :return:
        """
        assert isinstance(f1, PolyRingModP)
        assert isinstance(f2, PolyRingModP)
        assert isinstance(f3, PolyRingModP)

        return f1.fast_mul(f2) + f3

    def __mul__(self, other):
        """
            calculate the multiplication of two poly
        Notation:
            we have ck = \sigma aibj (i+j=k mod N)
        :param other: other polynomial
        :return: multiplication result
        """
        # warnings.warn("if changes to convolution, faster or not?")

        if isinstance(other, PolyRingModP):
            assert self.N == other.N
            assert self.p == other.p

            c = np.zeros(self.N, dtype=DTYPE)
            for i in range(self.N):
                for j in range(self.N):
                    if i+j >= self.N:
                        c[i+j-self.N] += self.coefs[i] * other.coefs[j]
                    else:
                        c[i+j] += self.coefs[i] * other.coefs[j]

            return PolyRingModP(self.N, self.p, coefs=c)
        elif isinstance(other, int):
            return PolyRingModP(self.N, self.p, coefs=self.coefs*other)
        else:
            raise NotImplementedError

    def __rmul__(self, other):
        return self * other

    def __radd__(self, other):
        return self + other

    def __add__(self, other):
        if isinstance(other, int):
            c = np.copy(self.coefs)
            c[0] += other
            return PolyRingModP(self.N, self.p, coefs=c)
        elif isinstance(other, PolyRingModP):
            assert self.N == other.N
            assert self.p == other.p

            c = np.zeros(self.N, dtype=DTYPE)
            for i in range(self.N):
                c[i] = self.coefs[i] + other.coefs[i]
            return PolyRingModP(self.N, self.p, coefs=c)
        else:
            raise NotImplementedError

    def __neg__(self):
        return PolyRingModP(self.N, self.p, coefs=-self.coefs)

    def __sub__(self, other):
        return self + -other

    def __mod__(self, other):
        if isinstance(other, int):
            return PolyRingModP(self.N, self.p, coefs=self.coefs % other)
        else:
            raise NotImplementedError

    @staticmethod
    def div(p, coefs1, coefs2):
        # more detail in 6.3.4.1
        assert type(coefs1) is np.ndarray
        assert type(coefs2) is np.ndarray

        r = coefs1.copy()
        r = r % p
        # this seems slower than mod
        # r = r & (self.p-1)
        r_degree = PolyRingModP.get_degree(r)
        r = r[:r_degree+1]

        b = coefs2.copy()
        b = b % p
        # b = b & (other.p-1)
        b_degree = PolyRingModP.get_degree(b)
        b = b[:b_degree+1]

        u = PolyRingModP.exgcd_not_rec(b[-1], p)[0] % p
        q = np.zeros(max(r_degree-b_degree+1, 1), dtype=DTYPE)

        while r_degree >= b_degree:
            d = r_degree
            tmp = (1*u*r[d]) % p
            v = np.array([0]*(d-b_degree)+[tmp], dtype=DTYPE)
            # https://zhuanlan.zhihu.com/p/68073482
            # the convolution of two polynomials
            conv = np.convolve(v, b)
            r = (np.pad(r, (0, max(0, conv.size-r.size))) - np.pad(conv, (0, max(0, r.size-conv.size)))) % p
            q[d-b_degree] += tmp
            if not r_degree:
                break

            r_degree = PolyRingModP.get_degree(r)
            r = r[:r_degree + 1]


        return q % p, r

    @staticmethod
    def extended_euclidean_algorithm(p, coefs1, coefs2):
        """
            find u, v, d that satisfy au+bv=d and d = gcd(a,b)
        :param p: mod p
        :param a: polynomial a
        :param b: polynomial b
        :return: (u, v, d)
        """

        assert type(coefs1) is np.ndarray
        assert type(coefs2) is np.ndarray

        a = coefs1.copy()
        b = coefs2.copy()

        if len(b) == 1 and b[0] == 0:
            return 1, 0, a

        u = np.ones(1, dtype=DTYPE)
        d = a.copy()
        v1 = np.zeros(1, dtype=DTYPE)
        v3 = b.copy()

        while len(v3) > 1 or v3[0] != 0:
            q, t3 = PolyRingModP.div(p, d, v3)
            conv = np.convolve(q, v1)
            t1 = (np.pad(u, (0, max(0, conv.size-u.size))) - np.pad(conv, (0, max(0, u.size-conv.size)))) % p
            u = v1.copy()
            d = v3.copy()
            v1 = t1.copy()
            v3 = t3.copy()

        tmp = np.convolve(a, u)
        v, r = PolyRingModP.div(p, (np.pad(d, (0, max(0, tmp.size-d.size))) -
                                    np.pad(tmp, (0, max(0, d.size-tmp.size)))), b)
        assert len(r) == 1 and r[0] == 0

        u = PolyRingModP.remove_right_zeros(u)
        # print(u, v, d)
        return u, v, d

    @staticmethod
    def inverse_mod_p(N, p, a):
        """
            details can be found in 6.3.4.3
        :param p: the modular
        :param N: the degree of the polynomial
        :param a: the poly
        :return: polynomial b satisfy ab = 1 if b exists else return false
        """
        if type(a) is np.ndarray:
            b = np.zeros(N+1, dtype=DTYPE)
            b[-1] = 1
            b[0] = p - 1
            u, v, d = PolyRingModP.extended_euclidean_algorithm(p, a, b)
            if len(d) == 1:
                return (u * PolyRingModP.exgcd_not_rec(d[0], p)[0]) % p
            else:
                return None
        else:
            raise NotImplementedError

    @staticmethod
    def inverse_mod_p_e(N, p, a, e):
        """
            details can be found in 6.3.4.3
        :param e: p^e
        :return:
        """
        b = PolyRingModP.inverse_mod_p(N, p, a)
        b = np.pad(b, (0, max(0, N-b.size)))
        m = e
        p_e = p ** e
        if b is not None:
            n = 2
            while m > 0:
                b = 2 * b - PolyRingModP.mul(
                    N, a, PolyRingModP.mul(
                        N, b, b
                    )
                )
                b = b % p_e
                m = m // 2
                n = n * 2
            return b

    @staticmethod
    def centralized(a, q):
        if type(a) is np.ndarray:
            mid = q // 2
            b = a.copy()
            for i in range(b.size):
                if b[i] > mid:
                    b[i] = b[i] - q
            return b
        else:
            raise NotImplementedError





if __name__ == '__main__':
    poly1 = np.array([1, 1, 1, 0, 1, 0, 1], dtype=DTYPE)
    poly2 = np.array([1, 1, 0, 1, 1, 0, 0, 0, 1], dtype=DTYPE)
    print(PolyRingModP.extended_euclidean_algorithm(2, poly1, poly2))

    poly1 = np.array([1, 1, 1, 0, 1, 0, 1], dtype=DTYPE)
    poly2 = np.array([1, 1, 1, 1, 1, 1, 0, 1], dtype=DTYPE)
    poly3 = np.array([1, 1, 0, 1, 1, 0, 0, 0, 1], dtype=DTYPE)
    poly4 = np.array([0, 0, 1, 0, 0, 1], dtype=DTYPE)
    print(np.convolve(poly1, poly2) + np.convolve(poly3, poly4))