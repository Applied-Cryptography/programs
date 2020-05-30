from utils import logger, reverse_int

from typing import List, Set, Tuple


def power_mod(n: int, a: int, p: int) -> List[int]:
    """返回 x^n ≡ a (mod p) 的所有解"""
    result_list = []

    for i in range(p):
        if i ** n % p == a:
            result_list.append(i)

    return result_list


class EllipticCurveModP:
    """域 p 下的椭圆曲线 y^2 = x^3 + ax + b"""

    def __init__(self, a, b, p):
        self.a = a
        self.b = b
        self.p = p
        self.points_set = self.cal_all_points(info=False)

    def cal_all_points(self, info=True) -> Set[Tuple[int, int]]:
        # 点对 (x, y)
        points_set = set()
        for x in range(self.p):
            a = (x ** 3 + self.a * x + self.b) % self.p
            y_result = power_mod(2, a, self.p)
            if len(y_result):
                for y in y_result:
                    points_set.add((x, y))
                if info:
                    logger.info(f"x={x}, y^2={a} (mod {self.p}),  y ≡ {', '.join(map(str, y_result))} (mod {self.p})")
            else:
                if info:
                    logger.info(f"x={x}, y^2={a} (mod {self.p}), 无解")

        return points_set

    def same_add(self, p: Tuple[int, int]) -> Tuple[int, int]:
        """计算 p+p"""
        assert p in self.points_set

        t = (3 * (p[0] ** 2) + self.a) * reverse_int(2*p[1], self.p) % self.p
        logger.info(f"λ = {t}")

        x = (t ** 2 - 2 * p[0]) % self.p
        y = (t * (p[0]-x) - p[1]) % self.p

        logger.info(f"x = {x}, y = {y}")

        return x, y

    def diff_add(self, p1: Tuple[int, int], p2: Tuple[int, int]) -> Tuple[int, int]:
        """计算 p1+p2 且 p1 != p2"""
        assert p1 in self.points_set
        assert p2 in self.points_set

        t = ((p2[1] - p1[1]) * reverse_int(p2[0]-p1[0], self.p)) % self.p
        logger.info(f"λ = {t}")

        x = (t ** 2 - p1[0] - p2[0]) % self.p
        y = (t * (p1[0]-x) - p1[1]) % self.p

        logger.info(f"x = {x}, y = {y}")

        return x, y


if __name__ == '__main__':
    test = EllipticCurveModP(1, 6, 11)
    test.cal_all_points()
    print(test.diff_add((7, 9), (8, 3)))
