from utils import is_prime, logger, reverse_int

from typing import Union, Tuple, Set, Optional


def judge_by_euler_formula(a: int, p: int) -> bool:
    """Determine whether a is the quadratic residue of module p by euler formula

    Args:
        a: an integer
        p: the module, must be prime

    Returns:
        True if a is the quadratic residue otherwise False
    """
    if a % p == 0:
        raise ValueError("a module p can not be 0")
    assert is_prime(p)

    if p == 2:
        return True

    if a == 1:
        return True
    if a == p - 1:
        return p % 4 == 1
    if a == 2:
        return not bool(((p ** 2 - 1) // 8) & 1)
    else:
        return pow(a, (p - 1) // 2, p) == 1


class SquareRootModP:
    """x^2 ≡ a  (mod p)的一般解法"""
    @staticmethod
    def find_non_residue(p: int) -> int:
        """找到模 p 的二次非剩余"""
        for i in range(2, p):
            if not judge_by_euler_formula(i, p):
                return i

        raise ValueError(f"{p} 没有二次非剩余")

    @staticmethod
    def divide_p_minus_one(p: int) -> Tuple[int, int]:
        """将 p-1 分解为 2^t * s, 其中 s 是奇数

        Returns:
            (t, s) where p-1 = 2^t * s
        """
        assert is_prime(p)

        t = 0
        s = p-1
        while s & 1 == 0:
            s //= 2
            t += 1

        return t, s

    @staticmethod
    def find_result(a: int, p: int) -> Optional[Set[int]]:
        """找出 x^2 ≡ a  (mod p) 的所有解"""
        assert is_prime(p)
        if pow(a, (p-1) // 2, p) != 1:
            logger.info(f"无解")
            return

        t, s = SquareRootModP.divide_p_minus_one(p)
        logger.info(f"{p-1} = 2^{t} * {s}")
        non_residue = SquareRootModP.find_non_residue(p)
        logger.info(f"{p} 的一个二次非剩余是 {non_residue}")
        b = pow(non_residue, s, p)
        logger.info(f"b = {non_residue}^{s} = {b} (mop {p})")

        result_set = set()

        x = pow(a, (s+1) // 2, p)
        reverse_a = reverse_int(a, p)
        t -= 1
        j = 0
        subscript_j = 0

        if t == 0:
            result_set.add(x)
            result_set.add(p - x)
            logger.info(f"{p} 是 4k+3 型素数，两个解为：{', '.join(map(str, result_set))}")
            return result_set

        logger.info(f"x_{t} = {a}^{(s+1) // 2} = {x}  (mod {p}),  a 的逆为 {reverse_a}")

        while t > 0:
            temp = pow(reverse_a * x**2, 2 ** (t-1), p)
            if temp % p == p-1:
                logger.info(f"(a^-1 * x^2_{t})^(2^{t-1}) ≡ -1 (mod {p})")
                j = 1
            else:
                logger.info(f"(a^-1 * x^2_{t})^2^{t-1} ≡ 1 (mod {p})")
                j = 0

            previous_x = x
            x = (x * b ** (j * 2 ** (subscript_j))) % p
            logger.info(f"j_{subscript_j} = {j}, x_{t-1} = x_{t} * b^(2^{subscript_j} * j_{subscript_j}) = {previous_x} * {b}^{2 ** subscript_j * j} = {x} (mod {p})")

            subscript_j += 1
            t -= 1

        result_set.add(x)
        result_set.add(p - x)
        logger.info(f"{p} 是 4k+1 型素数，两个解为：{', '.join(map(str, result_set))}")
        return result_set


def jacobi(n, m):
    """Determine weather n is not quadratic residue of m

    Args:
        n: an odd number
        m: the module

    Returns:
        -1 if n is not quadratic residue of m
        1 if can not decide
    """
    # todo
    # 需要先进行因式分解



if __name__ == '__main__':
    SquareRootModP.find_result(28, 37)
