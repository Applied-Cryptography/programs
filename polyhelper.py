from poly import PolyRingModP, DTYPE
from utils import logger

from typing import Iterable, List

import numpy as np


def print_poly(poly: List[int], extra_info: str = None, info=True) -> str:
    """按照指数从高到低打印多项式"""
    degree = len(poly) - 1

    poly_list = []
    for item in poly:
        if item == 0:
            degree -= 1
            continue
        item_add = f"{item}x^{degree}"
        if item == 1:
            item_add = f"x^{degree}"
        if degree == 0:
            item_add = item
        poly_list.append(f"{item_add}")
        degree -= 1

    if len(poly_list) == 0:
        poly_list = ["0"]

    poly_str = " + ".join(poly_list)

    if info:
        if extra_info:
            logger.info(f"{extra_info}: {poly_str}")
        else:
            logger.info(poly_str)

    return poly_str


def poly_mod(p: int, poly1: Iterable[int], poly2: Iterable[int]):
    """域 p 上多项式 poly1 模除 poly2

    Args:
        p: 模数
        poly1: 被除数，指数从高到低
        poly2: 除数，指数从高到低
    """

    p1 = np.array(poly1[::-1], dtype=DTYPE)
    p2 = np.array(poly2[::-1], dtype=DTYPE)

    q, r = PolyRingModP.div(p, p1, p2)

    print_poly(q[::-1], "不完全商")
    print_poly(r[::-1], "余数")

    return q[::-1], r[::-1]


def poly_mul(p: int, poly1: Iterable[int], poly2: Iterable[int], poly3: Iterable[int]):
    """域 p 下计算 poly1 * poly2 % poly3"""
    p1 = np.array(poly1[::-1], dtype=DTYPE)
    p2 = np.array(poly2[::-1], dtype=DTYPE)
    p3 = np.array(poly3[::-1], dtype=DTYPE)

    max_len = max(map(len, [p1, p2]))

    # 这里不需要 N
    p1 = PolyRingModP(N=max_len*3, p=p, coefs=p1)
    p2 = PolyRingModP(N=max_len*3, p=p, coefs=p2)

    mul = p1 * p2
    mul_result = print_poly(mul.coefs[::-1], "原乘积为")
    q, r = PolyRingModP.div(p, mul.coefs, p3)
    q_result = print_poly(q[::-1], info=False)
    r_result = print_poly(r[::-1], info=False)
    p3_result = print_poly(p3[::-1], info=False)

    logger.info(f"{mul_result} = ({q_result}) * ({p3_result}) + {r_result}")

    print_poly(r[::-1], "模之后为")

    return r[::-1]


def poly_reverse(p: int, poly1: Iterable[int], poly2: Iterable[int]):
    """域 p 下多项式 poly1 在模 poly2 下的逆"""

    p1 = np.array(poly1[::-1], dtype=DTYPE)
    p2 = np.array(poly2[::-1], dtype=DTYPE)

    s, t, d = PolyRingModP.extended_euclidean_algorithm(p, p1, p2)

    p1_result = print_poly(p1[::-1], info=False)
    p2_result = print_poly(p2[::-1], info=False)
    s_result = print_poly(s[::-1], info=False)
    t_result = print_poly(t[::-1], info=False)
    d_result = print_poly(d[::-1], info=False)

    logger.info("贝祖等式：")
    logger.info(
        f"({s_result}) * ({p1_result}) + ({t_result}) * ({p2_result}) = {d_result}"
    )
    logger.info("逆为：")
    logger.info(s_result)

    return s[::-1]


# todo
# 增加求所有的生成元


if __name__ == '__main__':

    a = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    b = [1, 0, 1, 1, 0, 0, 0]
    p = [1, 0, 0, 1, 1]

    # poly_reverse(2, a, p)
    poly_mod(2, a, p)