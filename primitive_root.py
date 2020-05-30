from utils import is_prime, factor, logger

from math import gcd
from typing import List


def find_one_primitive_root(p: int) -> int:
    """find one yuangen of prime p"""

    assert p >= 3
    assert is_prime(p)

    factor_dict = factor(p-1)
    # for logger
    temp = []
    for key, value in factor_dict.items():
        temp.append(f"{key}^{value}")

    logger.info(f"{p}-1 = {p-1} = {' * '.join(temp)}")

    factor_set = set(int((p-1)/k) for k in factor_dict.keys() )

    # todo
    # 这里可以将原根的证明表达出来
    for i in range(2, p):
        is_yuangen = True
        for j in factor_set:
            if i ** j % p == 1:
                is_yuangen = False
                break
        if is_yuangen:
            logger.info(f"{p} 的一个原根是 {i}")
            return i


def find_all_primitive_root(p: int) -> List[int]:
    first_primitive_root = find_one_primitive_root(p)

    primitive_root_list = [first_primitive_root]
    for i in range(2, p):
        if gcd(i, p-1) == 1:
            new_primitive_root = first_primitive_root ** i % p
            logger.info(f"{first_primitive_root}^{i} ≡ {new_primitive_root}  (mod {p})")
            primitive_root_list.append(new_primitive_root)

    primitive_root_list.sort()
    logger.info(f"p的所有原根是: {primitive_root_list}")

    return primitive_root_list


def print_zhibiaobiao(primitive_root: int, p: int) -> None:
    assert primitive_root in find_all_primitive_root(p)

    zhibiaobiao = []
    for i in range(0, p-1):
        zhibiaobiao.append((i, primitive_root ** i % p))

    zhibiaobiao.sort(key=lambda x: x[1])

    logger.info("指标表：")
    for i, j in zhibiaobiao:
        logger.info(f"{j} ≡ {primitive_root}^{i}  (mod {p})")


if __name__ == '__main__':
    # find_one_primitive_root(41)
    print_zhibiaobiao(6, 41)