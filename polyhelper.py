from poly import PolyRingModP, DTYPE
from utils import logger, factor

from typing import Iterable, List, Tuple, Set, Optional
from itertools import product
import math

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


def poly_mod(p: int, poly1: Iterable[int], poly2: Iterable[int], info=True):
    """域 p 上多项式 poly1 模除 poly2

    Args:
        p: 模数
        poly1: 被除数，指数从高到低
        poly2: 除数，指数从高到低
    """

    p1 = np.array(poly1[::-1], dtype=DTYPE)
    p2 = np.array(poly2[::-1], dtype=DTYPE)

    q, r = PolyRingModP.div(p, p1, p2)

    q_str = print_poly(q[::-1], info=False)
    r_str = print_poly(r[::-1], info=False)
    p1_str = print_poly(p1[::-1], info=False)
    p2_str = print_poly(p2[::-1], info=False)

    if info:
        logger.info(f"{p1_str} ≡ ({p2_str}) * ({q_str}) + {r_str}")

    return q[::-1], r[::-1]


def poly_add(p: int, poly1: Iterable[int], poly2: Iterable[int], poly3: Iterable[int], info=True):
    """域 p 下计算 poly1 * poly2 % poly3"""
    p1 = np.array(poly1[::-1], dtype=DTYPE)
    p2 = np.array(poly2[::-1], dtype=DTYPE)
    p3 = np.array(poly3[::-1], dtype=DTYPE)

    max_len = max(map(len, [p1, p2]))

    # 这里不需要 N
    p1 = PolyRingModP(N=max_len*3, p=p, coefs=p1)
    p2 = PolyRingModP(N=max_len*3, p=p, coefs=p2)

    add = p1 + p2
    add_result = print_poly(add.coefs[::-1], "原和为", info=info)
    q, r = PolyRingModP.div(p, add.coefs, p3)
    q_result = print_poly(q[::-1], info=False)
    r_result = print_poly(r[::-1], info=False)
    p3_result = print_poly(p3[::-1], info=False)

    if info:
        logger.info(f"{add_result} = ({q_result}) * ({p3_result}) + {r_result}")

    print_poly(r[::-1], "模之后为", info=info)

    return r[::-1]


def poly_mul(p: int, poly1: Iterable[int], poly2: Iterable[int], poly3: Iterable[int], info=True):
    """域 p 下计算 poly1 * poly2 % poly3"""
    p1 = np.array(poly1[::-1], dtype=DTYPE)
    p2 = np.array(poly2[::-1], dtype=DTYPE)
    p3 = np.array(poly3[::-1], dtype=DTYPE)

    max_len = max(map(len, [p1, p2]))

    # 这里不需要 N
    p1 = PolyRingModP(N=max_len*3, p=p, coefs=p1)
    p2 = PolyRingModP(N=max_len*3, p=p, coefs=p2)

    mul = p1 * p2
    mul_result = print_poly(mul.coefs[::-1], "原乘积为", info=info)
    q, r = PolyRingModP.div(p, mul.coefs, p3)
    q_result = print_poly(q[::-1], info=False)
    r_result = print_poly(r[::-1], info=False)
    p3_result = print_poly(p3[::-1], info=False)

    if info:
        logger.info(f"{mul_result} = ({q_result}) * ({p3_result}) + {r_result}")

    print_poly(r[::-1], "模之后为", info=info)

    return r[::-1]


def poly_modular_exponentiation(p: int, poly1: List[int], e: int, poly2: List[int]):
    """域 p 下计算 poly1^e % poly2"""
    assert e >= 1

    p1 = np.array(poly1[::-1], dtype=DTYPE)

    max_len = max((len(poly1)-1)*e+1, len(poly2))

    p1 = PolyRingModP(N=max_len*3, p=p, coefs=p1)
    tmp_p = PolyRingModP(N=max_len*3, p=p, coefs=p1.coefs.copy())

    for i in range(e-1):
        tmp_p = tmp_p * p1
        tmp_p = tmp_p % p

    return poly_mod(p, tmp_p.coefs[::-1], poly2)


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


class IrreduciblePolyModp:
    """GF(p) 下不可约多项式的判断"""
    def __init__(self, p=2):
        # 两个初始的不可约多项式 x, x+1
        self.p = p
        self.irreducible_dict = {
            # 1: [(1, 0), (1, 1)]
        }

    @staticmethod
    def _generate_all_n_degree_poly(p: int, n: int) -> Iterable:
        # 有空的话可以试试自己实现(回溯法应该可以)
        for poly in product(range(p), repeat=n):
            yield 1, *poly

    @staticmethod
    def _check_zero_order_poly(poly: Set[int]) -> bool:
        degree = len(poly)

        for coef in poly:
            if coef != 0:
                break
            degree -= 1

        return degree <= 1

    def generate_irreducible_poly(self, n, info=True) -> None:
        """产生所有的 n 阶不可约多项式"""
        degree = 1
        while degree <= n:
            for poly in self._generate_all_n_degree_poly(self.p, degree):
                if self._check_zero_order_poly(poly):
                    continue
                if self.irreducible_dict.get(degree) and poly in self.irreducible_dict[degree]:
                    continue

                is_irreducible = True
                for key, value in self.irreducible_dict.items():
                    if degree != 1 and key * 2 > degree:
                        break

                    for irreducible_poly in value:
                        r = poly_mod(self.p, poly, irreducible_poly, info=False)[1]
                        if len(r) == 1 and r == 0:
                            is_irreducible = False

                    if not is_irreducible:
                        break
                if is_irreducible:
                    if info:
                        print_poly(poly, "不可约多项式")
                    try:
                        self.irreducible_dict[degree].add(poly)
                    except KeyError:
                        self.irreducible_dict[degree] = {poly}

            degree += 1

    def is_irreducible_poly(self, poly: Tuple) -> bool:
        assert len(poly) >= 2

        self.generate_irreducible_poly(len(poly) // 2)

        for _, value in self.irreducible_dict.items():
            for irreducible_poly in value:
                r = poly_mod(self.p, poly, irreducible_poly)[1]
                if len(r) == 1 and r == 0:
                    logger.info("该多项式可约")
                    return False

        logger.info("该多项式不可约")
        return True


def is_generator(p: int, poly1: List[int], poly2: List[int]) -> bool:
    """判断 poly1 是否是扩域 GF(p)/poly2 的生成元"""
    assert IrreduciblePolyModp(p).is_irreducible_poly(tuple(poly2))

    order = p ** (len(poly2)-1)
    prime_factors = sorted(factor(order-1).keys())
    logger.info(f"{order-1} 的素因子为：{', '.join(map(str, prime_factors))}")
    exp_list = [(order-1) // prime_factor for prime_factor in prime_factors]

    for exp in exp_list:
        _, r = poly_modular_exponentiation(p, poly1, exp, poly2)
        if len(r) == 1 and r[0] == 1:
            logger.info(f"{print_poly(poly1, info=False)} 不是生成元")
            return False

    logger.info(f"{print_poly(poly1, info=False)} 是生成元")
    return True


def get_all_generator(p: int, poly1: List[int], poly2: List[int]):
    # ploy1其中一个生成元，p不可约多项式
    assert is_generator(p, poly1, poly2)
    num = len(poly2)-1
    num = p ** num - 1
    for i in range(1, num):
        if math.gcd(i, num) == 1:
            poly_modular_exponentiation(p, poly1, i, poly2)


# 求所有的n阶以下所有的不可约多项式
def get_all_irreducible_poly(p: int, n: int):
    IrreduciblePolyModp(p).generate_irreducible_poly(n)


def construct_gf_n(n: int, k_irreducible_poly: List[int] = None, int_mode=True) -> None:
    """构造 GF(n)"""
    from prettytable import PrettyTable
    pt = PrettyTable()

    def poly_to_integer(p: int, _poly: List[int]):
        degree = len(_poly) - 1
        to_integer: int = 0
        for j in range(degree+1):
            to_integer += p ** (degree-j) * _poly[j]

        return to_integer

    factor_dict = factor(n)
    factor_keys = list(factor_dict.keys())
    assert len(factor_keys) == 1

    p = factor_keys[0]
    k = factor_dict[p]

    if not k_irreducible_poly:
        # 找一个 k 次不可约多项式
        irreducible_poly_helper = IrreduciblePolyModp(p)
        irreducible_poly_helper.generate_irreducible_poly(k, info=False)
        k_irreducible_poly = sorted(list(irreducible_poly_helper.irreducible_dict[k]))[0]
        # 清空原来的不可约多项式集合
        irreducible_poly_helper.irreducible_dict = dict()
        irreducible_poly_helper.is_irreducible_poly(k_irreducible_poly)
    else:
        assert len(k_irreducible_poly) == k+1

    print_poly(k_irreducible_poly, f"{k} 次不可约多项式:")

    # 表头元素
    table_header = []
    for poly in product(range(p), repeat=k):
        table_header.append(poly)
    if int_mode:
        columns = [""] + (list(map(lambda x: str(poly_to_integer(p, x)), table_header)))
    else:
        columns = [""] + list(map(lambda x: print_poly(x, info=False), table_header))
    pt.field_names = columns

    logger.info("=" * 100)
    logger.info("加法表为：")
    for i, poly in enumerate(table_header):
        row = [poly]
        for addend_poly in table_header:
            row.append(poly_add(p, poly, addend_poly, k_irreducible_poly, info=False))
        if int_mode:
            pt.add_row(list(map(lambda x: poly_to_integer(p, x), row)))
        else:
            pt.add_row(list(map(lambda x: print_poly(x, info=False), row)))
    logger.info(pt)

    # 清除行
    for i in range(len(table_header)-1, -1, -1):
        pt.del_row(i)

    # 复制一遍算了
    logger.info("=" * 100)
    logger.info("乘法表为：")
    for i, poly in enumerate(table_header):
        row = [poly]
        for multiplicand_poly in table_header:
            row.append(poly_mul(p, poly, multiplicand_poly, k_irreducible_poly, info=False))
        if int_mode:
            pt.add_row(list(map(lambda x: poly_to_integer(p, x), row)))
        else:
            pt.add_row(list(map(lambda x: print_poly(x, info=False), row)))

    logger.info(pt)




if __name__ == '__main__':
    construct_gf_n(8)

    construct_gf_n(4, [1, 1, 1], int_mode=False)
