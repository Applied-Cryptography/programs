"""
below are some functions implement simple "RSA" and "exponentiation by squaring"
@author lyl
@time 2020/3/4
"""

from typing import List, Iterable
from random import randint

from utils import gen_nbit_prime


class RSA:
    """RSA example

    Description:
        this implement 2.4.5 in chapter2
    """

    def __init__(self, nbit, e=65537):

        p = gen_nbit_prime(nbit)
        q = gen_nbit_prime(nbit)
        n = p * q
        phi_n = (p-1) * (q-1)
        d, *_ = self.extended_euclidean(e, phi_n)
        d %= phi_n
        self.private_key = {
            "p": p,
            "q": q,
            "phi_n": phi_n,
            "d": d
        }

        self.public_key = {
            "n": n,
            "e": e
        }

    def encrypt(self, a):
        """encrypt plaintext a with public key e"""
        e = self.public_key.get("e")
        n = self.public_key.get("n")
        return pow(a, e, n)

    def decrypt(self, c):
        d = self.private_key.get("d")
        n = self.public_key.get("n")
        return pow(c, d, n)

    @staticmethod
    def extended_euclidean(a, b):
        """use extended euclidean algorithm to calculate s, t, k, s.t. sa + tb = k = (a, b)"""

        quotients = []
        q, r1 = divmod(a, b)
        r2 = b
        quotients.append(q)
        while r1 > 0:
            q, r = divmod(r2, r1)
            r2 = r1
            r1 = r
            quotients.append(q)

        s: list = [1, 0]
        t: list = [0, 1]
        for i in range(len(quotients) - 1):
            s.append(quotients[i] * -s[i+1] + s[i])
            t.append(quotients[i] * -t[i+1] + t[i])

        return s[-1], t[-1], r2


if __name__ == '__main__':
    # print(exponentiation_by_squaring(137, 113, 187))
    #
    rsa = RSA(512)
    plaintext = 1234567890
    encrypt = rsa.encrypt(plaintext)
    decrypt = rsa.decrypt(encrypt)
    print(f"plaintext: {plaintext}")
    print(f"ciphertext: {encrypt}")
    print(f"decryptext: {decrypt}")

    # prime_engine = GeneratePrime()
    # print(prime_engine.is_prime(prime_engine.gen_nbit_prime(1024) * prime_engine.gen_nbit_prime(1024))))