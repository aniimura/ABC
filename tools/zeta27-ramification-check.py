# -*- coding: utf-8 -*-
"""厳密整数演算による ℚ₃(ζ₂₇)/ℚ₃ の分岐データの検算。

★丸めを一切使わない。すべて Python の任意精度整数。

設定:
  M = ℚ₃(ζ₂₇)、π = ζ₂₇ − 1、E(π) = Φ₂₇(1+π) は π について Eisenstein、deg 18。
  𝒪_M = ℤ₃[π]、{π^0,…,π^17} が整基底。
  ⇒ y = Σ b_i π^i (b_i ∈ ℤ) について v_M(y) = min_i (18·v₃(b_i) + i)。
     （18·v₃(b_i) + i は i mod 18 で相異なるので min が実現する。）

  σ_a : ζ₂₇ ↦ ζ₂₇^a、すなわち π ↦ (1+π)^a − 1  (a ∈ (ℤ/27)^×)。
  i_G(σ_a) = v_M(σ_a(π) − π)。

検算対象（WildDescentDistanceOnly.lean:60-65 の字面）:
  e_M = 18
  G_0 = G, G_1 = G_2 = H, G_3..G_8 = Gal(M/ℚ₃(ζ₉)), G_9 = 1
  位数 9 の σ は i_G(σ) = 3、位数 3 の τ は i_G(τ) = 9、層 M/ℚ₃(ζ₉) の跳び i = 8
"""

N = 18  # deg E = φ(27)

# Φ₂₇(t) = t^18 + t^9 + 1
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


def shift(p, k):
    return [0] * k + list(p)


def compose_shift1(coeffs):
    """f(t) の係数から f(1+x) の係数を返す（整数、厳密）。"""
    res = [0]
    # Horner: f(1+x) = ((a_n (1+x) + a_{n-1})(1+x) + ...)
    for a in reversed(coeffs):
        res = poly_mul(res, [1, 1])
        res = poly_add(res, [a])
    return res


E = compose_shift1(PHI27)  # E(x) = Φ₂₇(1+x)
while len(E) > 1 and E[-1] == 0:
    E.pop()
assert len(E) == N + 1, (len(E), E)
assert E[-1] == 1, E[-1]


def v3(n):
    if n == 0:
        return None
    k = 0
    while n % 3 == 0:
        n //= 3
        k += 1
    return k


def reduce_mod_E(p):
    """deg < N に落とす（E は monic なので整数のまま）。"""
    p = list(p)
    for d in range(len(p) - 1, N - 1, -1):
        c = p[d]
        if c == 0:
            continue
        p[d] = 0
        for j in range(N):
            p[d - N + j] -= c * E[j]
    p = p[:N] + [0] * max(0, N - len(p))
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
    return best  # None = 0


PI = [0] * N
PI[1] = 1
ONE = [0] * N
ONE[0] = 1


def sigma_pi(a):
    """σ_a(π) = (1+π)^a − 1。"""
    base = poly_add(ONE, PI)
    acc = list(ONE)
    for _ in range(a):
        acc = mul(acc, base)
    return sub(acc, ONE)


def main():
    units = [a for a in range(1, 27) if a % 3 != 0]
    print("E(x) = Phi27(1+x) coefficients (low->high):", E)
    print("Eisenstein check: v3(E[0]) =", v3(E[0]),
          " all v3(E[i])>=1 for 0<i<18:",
          all(E[i] % 3 == 0 for i in range(1, N)))
    print()
    rows = []
    for a in units:
        sp = sigma_pi(a)
        i_sigma = vM(sub(sp, PI))
        # order of a in (Z/27)^*
        o = 1
        x = a % 27
        while x != 1:
            x = (x * a) % 27
            o += 1
        rows.append((a, o, i_sigma))
    print("a  ord(a in (Z/27)^*)  i_G(sigma_a) = v_M(sigma_a(pi) - pi)")
    for a, o, i in rows:
        print(f"{a:3d}  {o:3d}   {str(i)}")
    print()
    # ramification groups G_n = { sigma : i(sigma) >= n+1 }
    print("n   |G_n|   members(a)")
    for n in range(0, 12):
        mem = [a for a, o, i in rows if (i is not None and i >= n + 1)]
        mem_with_id = [1] + [a for a in mem if a != 1]
        print(f"{n:2d}  {len(mem_with_id):3d}   {sorted(set(mem_with_id))}")
    print()
    # Gal(M/K) with K = Q3(zeta3): a = 1 mod 3
    HK = [a for a in units if a % 3 == 1]
    print("Gal(M/K), K = Q3(zeta3):  a = 1 mod 3 :", sorted(HK), " size =", len(HK))
    # Gal(M/F) with F = Q3(zeta9): a = 1 mod 9
    HF = [a for a in units if a % 9 == 1]
    print("Gal(M/F), F = Q3(zeta9):  a = 1 mod 9 :", sorted(HF), " size =", len(HF))
    print()
    d = {a: i for a, o, i in rows}
    print("i_G on Gal(M/K):", {a: d[a] for a in sorted(HK)})
    print("i_G on Gal(M/F):", {a: d[a] for a in sorted(HF)})


if __name__ == "__main__":
    main()
