# -*- coding: utf-8 -*-
"""厳密整数演算による ℚ₃(ζ₂₇) の「実在する x」の再導出。

★WildDescentDistanceOnly.lean:57 は「スクリプトは報告に添付」と書いており、
★x もスクリプトもリポジトリに無い。本スクリプトはそれを再導出する。

設定は tools/zeta27-ramification-check.py と同じ（π = ζ₂₇ − 1、𝒪_M = ℤ₃[π]）。

検算対象（WildDescentDistanceOnly.lean:63-66）:
  x = 素元 のとき   v(σx−x) = 3,  v(τx−x) = 9
  ★実在する x では  v(σx−x) = 11, v(τx−x) = 15,  d(x,ℚ₃(ζ₉)) は v_M で 7
  σ = 位数 9 の元、τ = 位数 3 の元。
"""

import random
from fractions import Fraction

N = 18

PHI27 = [0] * 19
PHI27[18] = 1
PHI27[9] = 1
PHI27[0] = 1


def poly_mul(p, q):
    r = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a == 0:
            continue
        for j, b in enumerate(q):
            if b == 0:
                continue
            r[i + j] += a * b
    return r


def poly_add(p, q):
    n = max(len(p), len(q))
    r = [0] * n
    for i, a in enumerate(p):
        r[i] += a
    for i, b in enumerate(q):
        r[i] += b
    return r


def compose_shift1(coeffs):
    res = [0]
    for a in reversed(coeffs):
        res = poly_mul(res, [1, 1])
        res = poly_add(res, [a])
    return res


E = compose_shift1(PHI27)
while len(E) > 1 and E[-1] == 0:
    E.pop()
assert len(E) == N + 1 and E[-1] == 1


def v3(n):
    if isinstance(n, Fraction):
        if n == 0:
            return None
        return v3int(n.numerator) - v3int(n.denominator)
    return v3int(n)


def v3int(n):
    if n == 0:
        return None
    n = abs(n)
    k = 0
    while n % 3 == 0:
        n //= 3
        k += 1
    return k


def reduce_mod_E(p):
    p = list(p)
    for d in range(len(p) - 1, N - 1, -1):
        c = p[d]
        if c == 0:
            continue
        p[d] = 0
        for j in range(N):
            p[d - N + j] -= c * E[j]
    p = p + [0] * max(0, N - len(p))
    return p[:N]


def mul(a, b):
    return reduce_mod_E(poly_mul(a, b))


def sub(a, b):
    return [x - y for x, y in zip(a, b)]


def vM(y):
    best = None
    for i, b in enumerate(y):
        if b == 0:
            continue
        w = 18 * v3(b) + i
        if best is None or w < best:
            best = w
    return best


PI = [0] * N
PI[1] = 1
ONE = [0] * N
ONE[0] = 1


def pow_mod(base, k):
    acc = list(ONE)
    for _ in range(k):
        acc = mul(acc, base)
    return acc


def sigma_pi(a):
    return sub(pow_mod(poly_add(ONE, PI)[:N] if len(poly_add(ONE, PI)) <= N
                       else reduce_mod_E(poly_add(ONE, PI)), a), ONE)


def apply_sigma(a, y):
    """σ_a(y) —— y = Σ b_j π^j に σ_a(π) を代入。"""
    sp = sigma_pi(a)
    acc = [0] * N
    term = list(ONE)
    for j in range(N):
        if y[j] != 0:
            acc = [u + y[j] * v for u, v in zip(acc, term)]
        term = mul(term, sp)
    return acc


# ---- 𝒪_M = 𝒪_F[π] の基底 {μ^s π^i}, μ = ζ₉ − 1 = (1+π)^9 − 1, s<6, i<3 ----

MU = sub(pow_mod(reduce_mod_E(poly_add(ONE, PI)), 3), ONE)   # zeta9 - 1, v_M = 3
MUK = sub(pow_mod(reduce_mod_E(poly_add(ONE, PI)), 9), ONE)  # zeta3 - 1, v_M = 9


def basis_F_matrix():
    cols = []
    for s in range(6):
        ms = pow_mod(MU, s)
        for i in range(3):
            cols.append(mul(ms, pow_mod(PI, i)))
    return cols  # 18 個、各々長さ 18


COLS = basis_F_matrix()


def solve_coeffs(y):
    """y を {μ^s π^i} で展開する（厳密有理数）。返り値 c[s][i]。"""
    n = 18
    A = [[Fraction(COLS[k][r]) for k in range(n)] for r in range(n)]
    b = [Fraction(y[r]) for r in range(n)]
    # Gauss
    for col in range(n):
        piv = None
        for r in range(col, n):
            if A[r][col] != 0:
                piv = r
                break
        assert piv is not None, "singular"
        A[col], A[piv] = A[piv], A[col]
        b[col], b[piv] = b[piv], b[col]
        pv = A[col][col]
        A[col] = [v / pv for v in A[col]]
        b[col] = b[col] / pv
        for r in range(n):
            if r != col and A[r][col] != 0:
                f = A[r][col]
                A[r] = [u - f * v for u, v in zip(A[r], A[col])]
                b[r] = b[r] - f * b[col]
    c = [[b[s * 3 + i] for i in range(3)] for s in range(6)]
    return c


def vM_of_F_part(cs):
    """f = Σ_s c_s μ^s ∈ 𝒪_F について v_M(f) = min_s (18 v3(c_s) + 3 s)。"""
    best = None
    for s, cc in enumerate(cs):
        if cc == 0:
            continue
        w = 18 * v3(cc) + 3 * s
        if best is None or w < best:
            best = w
    return best


def basis_K_matrix():
    cols = []
    for s in range(2):
        ms = pow_mod(MUK, s)
        for i in range(9):
            cols.append(mul(ms, pow_mod(PI, i)))
    return cols


COLSK = basis_K_matrix()


def solve_gen(cols, y, ns, ni):
    n = 18
    A = [[Fraction(cols[k][r]) for k in range(n)] for r in range(n)]
    b = [Fraction(y[r]) for r in range(n)]
    for col in range(n):
        piv = None
        for r in range(col, n):
            if A[r][col] != 0:
                piv = r
                break
        assert piv is not None, "singular"
        A[col], A[piv] = A[piv], A[col]
        b[col], b[piv] = b[piv], b[col]
        pv = A[col][col]
        A[col] = [v / pv for v in A[col]]
        b[col] = b[col] / pv
        for r in range(n):
            if r != col and A[r][col] != 0:
                f = A[r][col]
                A[r] = [u - f * v for u, v in zip(A[r], A[col])]
                b[r] = b[r] - f * b[col]
    return [[b[s * ni + i] for i in range(ni)] for s in range(ns)]


def dist_to_K(y):
    c = solve_gen(COLSK, y, 2, 9)
    out = []
    for i in range(1, 9):
        best = None
        for s in range(2):
            cc = c[s][i]
            if cc == 0:
                continue
            w = 18 * v3(cc) + 9 * s
            if best is None or w < best:
                best = w
        if best is not None:
            out.append(best + i)
    return min(out) if out else None


def dist_to_F(y):
    """d(y, F) を v_M で。x = f_0 + f_1 π + f_2 π² ⇒ min_{i=1,2}(v_M(f_i)+i)。"""
    c = solve_coeffs(y)
    out = []
    for i in (1, 2):
        cs = [c[s][i] for s in range(6)]
        v = vM_of_F_part(cs)
        if v is not None:
            out.append(v + i)
    return min(out) if out else None


def main():
    random.seed(20260909)
    print("i(sigma_4) =", vM(sub(sigma_pi(4), PI)),
          "  i(sigma_10) =", vM(sub(sigma_pi(10), PI)))
    print("uniformizer x = pi :  v(sigma x - x) =",
          vM(sub(apply_sigma(4, PI), PI)),
          " v(tau x - x) =", vM(sub(apply_sigma(10, PI), PI)))
    print("  d(pi, F) =", dist_to_F(PI), "  d(pi, K) =", dist_to_K(PI))
    print()
    hits = {}
    trials = 20000
    for _ in range(trials):
        y = [random.randrange(0, 81) for _ in range(N)]
        vs = vM(sub(apply_sigma(4, y), y))
        vt = vM(sub(apply_sigma(10, y), y))
        key = (vs, vt)
        hits[key] = hits.get(key, 0) + 1
        if key == (11, 15) and "sample" not in hits:
            hits["sample"] = list(y)
    print(f"trials = {trials}   (v(sigma x - x), v(tau x - x)) の分布 (上位 12):")
    items = [(k, v) for k, v in hits.items() if isinstance(k, tuple)]
    items.sort(key=lambda kv: -kv[1])
    for k, v in items[:12]:
        print("   ", k, v)
    if "sample" in hits:
        y = hits["sample"]
        print()
        print("★(11,15) の実例 x (pi-basis, low->high):", y)
        print("   v(sigma x - x) =", vM(sub(apply_sigma(4, y), y)))
        print("   v(tau   x - x) =", vM(sub(apply_sigma(10, y), y)))
        print("   d(x, F=Q3(zeta9)) =", dist_to_F(y))
        print("   d(x, K=Q3(zeta3)) =", dist_to_K(y))
    else:
        print()
        print("(11,15) は出なかった")


if __name__ == "__main__":
    main()
