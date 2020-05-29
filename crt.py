"""
  包括中国剩余定理等一些函数
"""


import math
from typing import List

from utils import extended_euclidean, logger

from numpy import poly1d


def crt(b_list: List[int], m_list: List[int]) -> int:
    assert len(b_list) == len(m_list)

    # 每一个 m_i 必须互相互质
    for i in range(len(m_list)):
        for j in range(i+1, len(m_list)):
            if math.gcd(m_list[i], m_list[j]) != 1:
                raise ValueError("m_i 必须互相互质")

    m = math.prod(m_list)
    logger.info(f"M: {m}")
    upper_m_list = [m // i for i in m_list]
    logger.info(f"M_i: {upper_m_list}")
    upper_m_prime_list = [
        extended_euclidean(upper_m_list[i], m_list[i])[0] % m_list[i]
        for i in range(len(upper_m_list))
    ]
    logger.info(f"M_i': {upper_m_prime_list}")

    result = sum([
        b_list[i]*upper_m_list[i]*upper_m_prime_list[i]
        for i in range(len(b_list))
    ]) % m

    logger.info(f"result: {result}")

    return result


# fixme
# 应该试探所有的解
# 考虑不够周到
def poly_q_e(poly: List[int], q: int, e: int) -> int:
    """求 poly ≡ 0 模 p^e"""
    i = 0
    p = poly1d(poly)
    logger.info("poly:")
    logger.info(f"{p}")
    while i < q and (p(i) % q != 0):
        i += 1

    if i >= q:
        raise ValueError("无解")
    logger.info(f"x_1 = {i}")

    p_derivative = p.deriv()
    logger.info("poly derivative:")
    logger.info(f"{p_derivative}")

    p_derivative_i = p_derivative(i)
    p_derivative_i_reverse = extended_euclidean(p_derivative_i, q)[0] % q
    logger.info(f"f'(x_1)^-1 ≡ {p_derivative_i_reverse}  (mod {q})")
    if math.gcd(p_derivative_i_reverse, q) != 1:
        raise ValueError("无解")

    j = 1
    x_j = i
    while j < e:
        t_j = ((-p(x_j) // (q ** j)) * p_derivative_i_reverse) % q
        x_j = (x_j + t_j * (q ** j)) % (q ** (j+1))
        logger.info(f"t_{j+1} = {t_j}, x_{j+1} ≡ {x_j}  (mod {q**(j+1)})")

        j += 1

    return x_j



if __name__ == '__main__':
    # poly_q_e(
    #     [1, 5, 0, 9],
    #     3,
    #     3
    # )

    crt([0, 2, 2, 2, 2], [13, 3, 5 ,7 ,11])

