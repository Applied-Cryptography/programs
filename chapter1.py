import numpy as np
import matplotlib.pyplot as plt


class FindPrime:
    """compare eratosthenes sieve algorithm with normal method"""

    @staticmethod
    def eratosthenes(n):
        """use eratosthenes sieve to count all prime number"""
        if n <= 1:
            return 0
        sifter = np.ones(n+1, dtype=np.int8)
        sifter[0:2] = 0
        for i in range(2, int(np.ceil(np.sqrt(n)))+1):
            if sifter[i]:
                for j in range(2, n // i + 1):
                    sifter[i * j] = 0

        return sifter.sum()

    @staticmethod
    def normal(n):
        count = 0
        for i in range(2, n+1):
            if i == 2:
                count += 1
                continue
            flag = 1
            for j in range(2, int(np.ceil(np.sqrt(i)))+1):
                if i % j == 0:
                    flag = 0
                    break
            count += flag

        return count


if __name__ == '__main__':
    n = 10000
    x = np.arange(2, n)
    results = []
    for i in x:
        results.append(FindPrime.eratosthenes(i))
    y = x / np.log(x)
    plt.plot(x, results, label="prime")
    plt.plot(x, y, label="lnx")
    plt.legend(loc="upper left")
    plt.show()