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
    first_yuangen = find_one_primitive_root(p)

    yuangen_list = [first_yuangen]
    for i in range(2, p):
        if gcd(i, p-1) == 1:
            new_yuangen = first_yuangen ** i % p
            logger.info(f"{first_yuangen}^{i} ≡ {new_yuangen}  (mod {p})")
            yuangen_list.append(new_yuangen)

    yuangen_list.sort()
    logger.info(f"p的所有原根是: {yuangen_list}")

    return yuangen_list


def print_zhibiaobiao(yuangen: int, p: int) -> None:
    assert yuangen in find_all_primitive_root(p)

    zhibiaobiao = []
    for i in range(0, p-1):
        zhibiaobiao.append((i, yuangen ** i % p))

    zhibiaobiao.sort(key=lambda x: x[1])

    logger.info("指标表：")
    for i, j in zhibiaobiao:
        logger.info(f"{j} ≡ {yuangen}^{i}  (mod {p})")


if __name__ == '__main__':
    find_all_primitive_root(7)