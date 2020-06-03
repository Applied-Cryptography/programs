from utils import is_prime, factor, logger, reverse_int

from math import gcd
from typing import List, Dict, Optional


def find_one_primitive_root(p: int, info=True) -> int:
    """find one primitive_root of prime p"""
    assert p >= 3
    assert is_prime(p)

    factor_dict = factor(p-1)
    # for logger
    temp = []
    for key, value in factor_dict.items():
        temp.append(f"{key}^{value}")

    if info:
        logger.info(f"{p}-1 = {p-1} = {' * '.join(temp)}")

    factor_set = set(int((p-1)/k) for k in factor_dict.keys() )

    # todo
    # 这里可以将原根的证明表达出来
    for i in range(2, p):
        is_primitive_root = True
        for j in factor_set:
            if i ** j % p == 1:
                is_primitive_root = False
                break
        if is_primitive_root:
            if info:
                logger.info(f"{p} 的一个原根是 {i}")
            return i


def find_all_primitive_root(p: int, info=True) -> List[int]:
    first_primitive_root = find_one_primitive_root(p)

    primitive_root_list = [first_primitive_root]
    for i in range(2, p):
        if gcd(i, p-1) == 1:
            new_primitive_root = first_primitive_root ** i % p
            if info:
                logger.info(f"{first_primitive_root}^{i} ≡ {new_primitive_root}  (mod {p})")
            primitive_root_list.append(new_primitive_root)

    primitive_root_list.sort()
    if info:
        logger.info(f"p的所有原根是: {primitive_root_list}")

    return primitive_root_list


def construct_ind_table(primitive_root: int, p: int, info=True) -> Dict[int, int]:
    """求模 p 意义下，一个原根的指标表

    Returns:
        返回的形式为 Dict[a, b], primitive_root^b ≡ a  (mod p)
    """
    assert primitive_root in find_all_primitive_root(p)

    ind_table = dict()
    for i in range(0, p-1):
        b = primitive_root ** i % p
        ind_table[b] = i

    if info:
        logger.info("指标表：")
    for i, j in sorted(ind_table.items(), key=lambda x: x[0]):
        if info:
            logger.info(f"{primitive_root}^{j} ≡ {i}  (mod {p})")

    return ind_table


def one_order_congruence(a: int, b: int, m: int, info=True) -> Optional[List[int]]:
    """同余式 ax ≡ b  (mod m) 的解"""
    gcd_am = gcd(a, m)
    if b % gcd_am != 0:
        if info:
            logger.info(f"({a}, {m}) = {gcd_am} 不能整除 {b}")
        return

    if info:
        logger.info(f"原方程等价于：")
    a_tmp, b_tmp, m_tmp = map(lambda x: x // gcd_am, [a, b, m])
    if info:
        logger.info(f"{a_tmp}x ≡ {b_tmp}  (mod {m_tmp})")
    reverse_a_tmp = reverse_int(a_tmp, m_tmp)
    x_tmp = b_tmp * reverse_a_tmp % m_tmp
    if info:
        logger.info(f"该方程的解为：x ≡ {x_tmp} (mod {m_tmp})")

    result_list = []
    for i in range(gcd_am):
        result = x_tmp + i * m_tmp
        result_list.append(result)
    if info:
        logger.info(f"原方程的所有解为：{', '.join(map(str, result_list))}")

    return result_list


def higher_order_congruence(a: int, n: int, b: int, p: int) -> None:
    """求解 ax^n ≡ b  (mod p)"""
    assert is_prime(p)

    reverse_a = reverse_int(a, p)
    logger.info(f"{a} 在模 {p} 下的逆为 {reverse_a}")
    b = b * reverse_a % p

    primitive_root = find_one_primitive_root(p, info=False)
    logger.info(f"{p} 的一个原根为 {primitive_root}")

    ind_table = construct_ind_table(primitive_root, p, info=False)
    logger.info(f"ind {b} = {ind_table.get(b)}")

    ind_x_results = one_order_congruence(n, ind_table.get(b), p-1, info=False)
    if len(ind_x_results):
        logger.info(f"ind x的所有解为：{', '.join(map(str, ind_x_results))}")
    else:
        logger.info(f"ind x 无解！")

    logger.info(f"x 的所有解为：{', '.join(map(lambda x: str(primitive_root ** x % p), ind_x_results))}")

if __name__ == '__main__':
    # find_one_primitive_root(47)
    find_all_primitive_root(31)
    # construct_ind_table(5, 47)

    # higher_order_congruence(1, 5, 29, 47)