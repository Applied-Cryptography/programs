"""
some useful functions
"""


import logging
from typing import List, Dict
from random import randint
from collections import defaultdict
from math import gcd

logging.basicConfig(format="%(message)s", level=logging.INFO)
logger = logging.getLogger("logger")


def miller_robin_test(p: int, a: int) -> bool:
    """miller robin algorithm to test if p is prime

    More details can be found in wikipedia "https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test"
    """
    # fermat theorem: a^p = a (mod p) if p is prime
    # here we use builtin function for efficiency
    if pow(a, p, p) != a:
        return False

    n = p - 1
    # if x^2 = 1 (mod p) then, x = 1/p-1 (mod p)
    while n & 1:
        remainder = pow(a, n, p)
        if remainder not in (1, p-1):
            return False
        n //= 2

    return True


def is_prime(p) -> bool:
    if p & 1 == 0:
        return False
    if p in {2, 3, 5, 7, 11, 13, 17, 19, 23}:
        return True

    for judge_num in {2, 3, 5, 7, 11, 13, 17, 19, 23}:
        if not miller_robin_test(p, judge_num):
            return False

    return True


def gen_nbit_prime(n, threshold=1000):
    """generate nbit length prime"""

    # support n <= 1024
    assert n <= 1024

    lower_bound = 2 ** n
    upper_bound = 2 ** (n+1)

    num = 0
    while num < threshold:
        # keep odd
        p = randint(lower_bound, upper_bound) | 1
        if is_prime(p):
            return p
        num += 1
    print("fail to find prime")


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
        s.append(quotients[i] * -s[i + 1] + s[i])
        t.append(quotients[i] * -t[i + 1] + t[i])

    return s[-1], t[-1], r2


def reverse_int(a: int, m: int) -> int:
    """计算 a 在模 m 下的逆"""
    assert gcd(a, m) == 1

    return extended_euclidean(a, m)[0] % m


# todo
# 不详细
def exponentiation_by_squaring(b: int, e: int, m: int) -> int:
    """calculate non-negative integer k that satisfies b^e = k (mod m)

    Description:
        This function implement exponentiation mod in chapter2

    Args:
        b: the base integer
        e: the exponent integer
        m: the modulus

    Returns:
        the non-negative integer k that satisfies b^e = k (mod m) and 0 <= k < m
    """

    def integer_to_binary(s) -> List[int]:
        """change non-negative integer s to binary"""
        bin_array = []
        while s > 0:
            bin_array.append(s % 2)
            s //= 2

        # the binary is expressed as a_0 + a_1*2 + ... + a_n*2^n
        return bin_array

    bin_list = integer_to_binary(e)
    # previously calculate b, b^2, b^4, ... in the case of module m
    previous_b_list: list = [b % m]
    for i in range(len(bin_list) - 1):
        previous_b_list.append(previous_b_list[i] ** 2 % m)

    # calculate k
    k = b ** bin_list[0] % m
    for i in range(1, len(bin_list)):
        if bin_list[i] == 1:
            k = k * previous_b_list[i] % m
    # debug
    print(previous_b_list)
    print(bin_list)

    return k


def factor(n: int) -> Dict[int, int]:
    result = defaultdict(int)

    temp = n
    for i in range(2, int(-(-n ** 0.5 // 1)) + 1):
        if temp % i == 0:
            cnt = 0
            while temp % i == 0:
                cnt += 1
                temp //= i
            result[i] = cnt

    if temp != 1:
        result[temp] = 1

    if len(result) == 0:
        result[n] = 1

    return result



if __name__ == '__main__':
    exponentiation_by_squaring(31, 29, 113)

