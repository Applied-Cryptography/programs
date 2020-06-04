<h1 style="text-align: center;">
    常用函数汇总
</h1>

### 数论



#### 常用

- $(a, m) = 1$ 时，求模 $m$ 下 $a^{-1}$

  ```python
  # 位于 utils.py
  reverse_int(a: int, m: int) -> int
  ```

- $a^n \quad (mod \; m)$

  ```python
  # 位于 utils.py
  exponentiation_by_squaring(a: int, n: int, m: int) -> int
  ```

- $n$ 标准分解

  ```python
  # 位于 utils.py
  # 其中 key 表示质因数，value 表示该质因数的幂
  factor(n: int) -> Dict[int, int]
  ```

- $ax \equiv b \quad (mod \; m)$ 的解

  ```python
  # 位于 primitive_root.py
  one_order_congruence(a: int, b: int, m: int, info=True) -> Optional[List[int]]
  ```

  

- $a_i x \equiv b_i \quad (mod \; m_i)$ 在 $(m_i, m_j) = 1 \quad i \ne j$ 时的解

  ```python
  # 位于 crt.py
  crt(a_list: List[int], b_list: List[int], m_list: List[int]) -> int
  ```



#### 中国剩余定理

- 求多项式在模 $m$ 意义下的解：

  - 暴力枚举

  ```python
  # 位于 crt.py
  # 传入的 poly 为多项式按照x次数从高到底的系数
  poly_by_enumerate(poly: List[int], m: int, info=True) -> List[int]
  ```

  - 若 $m$ 为素数，可以尝试以下函数<strong style="color: red;">(不建议单独使用此函数，实现的不够全面)</strong>

  ```python
  # 位于 crt.py
  # 传入的 poly 为多项式按照x次数从高到底的系数
  poly_q_e(poly: List[int], q: int, e: int) -> int
  ```



#### 二次剩余

- $(\cfrac{a}{p})$ 

  ```python
  # 并未直接实现勒让德符号的计算
  # 可以用位于 square_root.py 的欧拉判别代替
  judge_by_euler_formula(a: int, p: int) -> bool
  ```

- $x^2 \equiv a \quad (mod \;  p)$ 的解

  ```python
  # 位于 square_root.py
  SquareRootModP.find_result(a: int, p: int) -> Optional[Set[int]]
  ```



#### 原根

- $p$ 的一个原根

  ```python
  # 位于 primitive_root.py
  find_one_primitive_root(p: int, info=True) -> int
  ```

- $p$ 的所有原根

  ```python
  # 位于 primitive_root.py
  find_all_primitive_root(p: int, info=True) -> List[int]
  ```

- 已知 $p$ 的一个原根 $primitive\_root$，计算指标表

  ```python
  # 位于 primitice_root.py
  # 返回的形式为 Dict[a, b], primitive_root^b ≡ a  (mod p)
  construct_ind_table(primitive_root: int, p: int, info=True) -> Dict[int, int]
  ```

- 高次同余式 $ax^n \equiv b \quad (mod \; p)$ 的解：

  ```python
  # 位于 primitive_root.py
  higher_order_congruence(a: int, n: int, b: int, p: int) -> None
  ```



### 群

​		所有和群中多项式相关的计算均位于 polyhelper.py

- 域 $p$ 上多项式 $poly1 \quad (mod \; poly2)$

  ```python
  poly_mod(p: int, poly1: Iterable[int], poly2: Iterable[int], info=True)
  ```

- 域 $p$ 下计算 $poly1 * poly2 \quad (mod \; poly3)$

  ```python
  poly_mul(p: int, poly1: Iterable[int], poly2: Iterable[int], poly3: Iterable[int])
  ```

- 域 $p$ 下计算 $poly1^e \quad (mod \; poly2)$

  ```python
  poly_modular_exponentiation(p: int, poly1: List[int], e: int, poly2: List[int])
  ```

- 域 $p$ 下多项式 $poly1$ 在模 $poly2$ 下的逆

  ```python
  poly_reverse(p: int, poly1: Iterable[int], poly2: Iterable[int])
  ```

- $GF(p)$ 下不可约多项式的判断

  ```python
  IrreduciblePolyModp(p).is_irreducible_poly(self, poly: Tuple) -> bool
  ```

- 判断 $poly1$ 是否是扩域 $F_p[x]/poly2$ 的生成元

  ```python
  is_generator(p: int, poly1: List[int], poly2: List[int]) -> bool
  ```



### 椭圆曲线

​		位于elliptic_curve.py

- 域 $p$ 下的椭圆曲线 $y^2 = x^3 + ax + b$

  ```python
  e = EllipticCurveModP(a, b, p)
  ```

  - 计算所有在椭圆上的点

    ```python
    e.call_all_points()
    ```

  - 计算 $p+p$ <strong style="color: red;">这里有点问题，应该只是横坐标相等</strong>

    ```python'
    e.same_add(self, p: Tuple[int, int]) -> Tuple[int, int]
    ```

  - 计算 $p_1 + p_2$，且 $p_1$ 的横坐标不等于 $p_2$ 的横坐标

    ```python
    e.diff_add(self, p1: Tuple[int, int], p2: Tuple[int, int]) -> Tuple[int, int]
    ```



### 离散数学

​		位于 truth_table.py

- 构建真值表

  ```python
  Solve("(pvq->q^r)->p^!r")
  ```

  

