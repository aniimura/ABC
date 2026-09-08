# -*- coding: utf-8 -*-
"""厳密整数演算による ℚ_p(ζ_{p^n})/ℚ_p の分岐群・跳びの検算（一般 p, n）。

★丸めを一切使わない。すべて Python の任意精度整数。

設定:
  L = ℚ_p(ζ_{p^n})、π = ζ_{p^n} − 1、E(π) = Φ_{p^n}(1+π) は Eisenstein、deg = φ(p^n)。
  𝒪_L = ℤ_p[π] が整基底 ⇒ v_L(Σ b_i π^i) = min_i (e·v_p(b_i) + i)、e = φ(p^n)。
  σ_a : π ↦ (1+π)^a − 1  (a ∈ (ℤ/p^n)^×)、i_G(σ_a) = v_L(σ_a(π) − π)。

出力:
  * 分岐群 G_u = { σ : i(σ) ≥ u+1 } の位数の列
  * 跳び（break）= G_u ≠ G_{u+1} となる u
  * 中間体 ℚ_p(ζ_{p^j}) に対応する部分群と、その層の跳び

検算対象:
  JumpFromValueGroup.lean:332 harith_zeta81 の u = (2, 8, 26)（p=3, n=4）
  GainedTowerDescent.lean:411 zeta27_sharp の L = 8
  GainedTowerDescent.lean:417 zeta81_sharp の L = 26
  GainedTowerDescent.lean:423 zeta16_sharp（p=2, n=4）
  WildDescentDistanceOnly.lean:60-62（p=3, n=3）—— tools/zeta27-ramification-check.py と同じ
"""


def cyclotomic_prime_power(p, n):
    """Φ_{p^n}(t) の係数（低次→高次）。Φ_{p^n}(t) = Σ_{k<p} t^{k·p^{n-1}}。"""
    q = p ** (n - 1)
    deg = q * (p - 1)
    c = [0] * (deg + 1)
    for k in range(p):
        c[k * q] = 1
    return c


def poly_mul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b):
            if y == 0:
                continue
            r[i + j] += x * y
    return r


def poly_add(a, b):
    n = max(len(a), len(b))
    r = [0] * n
    for i, x in enumerate(a):
        r[i] += x
    for i, y in enumerate(b):
        r[i] += y
    return r


def compose_shift1(coeffs):
    res = [0]
    for a in reversed(coeffs):
        res = poly_mul(res, [1, 1])
        res = poly_add(res, [a])
    return res


class Field:
    def __init__(self, p, n):
        self.p = p
        self.n = n
        self.q = p ** n
        phi = cyclotomic_prime_power(p, n)
        E = compose_shift1(phi)
        while len(E) > 1 and E[-1] == 0:
            E.pop()
        self.N = len(E) - 1           # e = φ(p^n)
        self.E = E
        assert E[-1] == 1
        # Eisenstein check
        assert self.vp(E[0]) == 1, ("E0 not Eisenstein", E[0])
        assert all(x % p == 0 for x in E[:self.N]), "middle coeffs"
        self.ONE = [0] * self.N
        self.ONE[0] = 1
        self.PI = [0] * self.N
        self.PI[1] = 1

    def vp(self, m):
        if m == 0:
            return None
        m = abs(m)
        k = 0
        while m % self.p == 0:
            m //= self.p
            k += 1
        return k

    def red(self, poly):
        N, E = self.N, self.E
        poly = list(poly)
        for d in range(len(poly) - 1, N - 1, -1):
            c = poly[d]
            if c == 0:
                continue
            poly[d] = 0
            for j in range(N):
                poly[d - N + j] -= c * E[j]
        poly = poly + [0] * max(0, N - len(poly))
        return poly[:N]

    def mul(self, a, b):
        return self.red(poly_mul(a, b))

    def sub(self, a, b):
        return [x - y for x, y in zip(a, b)]

    def v(self, y):
        best = None
        for i, b in enumerate(y):
            if b == 0:
                continue
            w = self.N * self.vp(b) + i
            if best is None or w < best:
                best = w
        return best

    def powm(self, base, k):
        acc = list(self.ONE)
        cur = list(base)
        while k > 0:
            if k & 1:
                acc = self.mul(acc, cur)
            cur = self.mul(cur, cur)
            k >>= 1
        return acc

    def sigma_pi(self, a):
        return self.sub(self.powm(poly_add(self.ONE, self.PI)[: self.N], a), self.ONE)

    def units(self):
        return [a for a in range(1, self.q) if a % self.p != 0]


def report(p, n):
    F = Field(p, n)
    print(f"=== L = Q_{p}(zeta_{p}^{n})   e = deg E = {F.N} ===")
    ivals = {}
    for a in F.units():
        ivals[a] = F.v(F.sub(F.sigma_pi(a), F.PI))
    sizes = []
    umax = (F.N // (p - 1)) * p + 4
    for u in range(0, umax):
        mem = [a for a in F.units() if a == 1 or (ivals[a] is not None and ivals[a] >= u + 1)]
        sizes.append((u, len(mem)))
    breaks = [u for (u, s), (u2, s2) in zip(sizes, sizes[1:]) if s != s2]
    print("  |G_u| (u = 0,1,2,...):", [s for _, s in sizes[: min(len(sizes), 30)]])
    print("  ★breaks (G_u != G_{u+1}) :", breaks)
    wild = [u for u in breaks if u >= 1]
    print("  ★wild breaks             :", wild)
    for j in range(1, n):
        sub = [a for a in F.units() if a % (p ** j) == 1]
        iv = [ivals[a] for a in sub if a != 1]
        if iv:
            print(f"  Gal(L / Q_{p}(zeta_{p}^{j})) : size {len(sub)},  "
                  f"min i = {min(iv)}  ⇒ 層の跳び = {min(iv) - 1}")
    print()


def main():
    report(3, 3)   # Q3(zeta27) —— WildDescentDistanceOnly.lean:60-62
    report(3, 4)   # Q3(zeta81) —— harith_zeta81 の u = (2,8,26)、zeta81_sharp の L = 26
    report(2, 4)   # Q2(zeta16) —— zeta16_sharp


if __name__ == "__main__":
    main()
