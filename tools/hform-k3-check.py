# -*- coding: utf-8 -*-
"""★`hform` の破れの「鋭い定数」を k = 2 と k = 3 の両方で測り、法則があるかを見る。

tools/hform-population-check.py は L = ℚ₃(ζ₂₇)（k = 2, e = 18）だけを測った:
  ★損失指数 eps − d(x,E₁) の最大は **4**（140,000 回の評価で 4 を超えず）。
  ★i₁ = 2 なので hform は 99.7% で成り立ち、稀に指数で 2 だけ破れる。

本スクリプトは同じ測定を **一般の (p, n)** で行い、k = 3（L = ℚ₃(ζ₈₁), e = 54）を測る。
  L = ℚ_p(ζ_{p^n}),  π = ζ_{p^n} − 1,  e = φ(p^n)
  E₁ = ℚ_p(ζ_{p^{n−1}}),  μ = (1+π)^p − 1 は E₁ の素元（v_L(μ) = p）
  𝒪_L = 𝒪_{E₁}[π] の ℤ_p 基底 {μ^s π^i}, s < e/p, i < p
  H_K = Gal(L/ℚ_p(ζ_p)) = {σ_a : a ≡ 1 (mod p)}
  eps(x)  = min_{σ ∈ H_K, σ≠1} v_L(σx − x)
  d(x,E₁) = min_{1≤i<p} ( v_L(f_i) + i ),   x = Σ_i f_i π^i
  ★損失指数 = eps − d(x,E₁)、hform は「損失指数 ≤ i₁」を要求する。

★k = 2 を再現できることを自己検査で確かめてから k = 3 を測る。
"""

import importlib.util
import random
import sys
import time
from fractions import Fraction

T = r"D:/Math_ABC3/tools/zeta-tower-check.py"
_spec = importlib.util.spec_from_file_location("ztow3", T)
zt = importlib.util.module_from_spec(_spec)
sys.modules["ztow3"] = zt
_spec.loader.exec_module(zt)


def _gcd(a, b):
    while b:
        a, b = b, a % b
    return a


class Layer:
    """L = ℚ_p(ζ_{p^n}) と、その一段下 E₁ = ℚ_p(ζ_{p^{n−1}}) への降下の測定器。"""

    def __init__(self, p, n):
        self.p, self.n = p, n
        self.F = zt.Field(p, n)
        F = self.F
        self.N = F.N
        self.mu = F.sub(F.powm(zt.poly_add(F.ONE, F.PI)[: F.N], p), F.ONE)
        assert F.v(self.mu) == p, ("v(mu)", F.v(self.mu))
        self.S = F.N // p                       # 𝒪_{E₁} の ℤ_p 階数
        self.HK = [a for a in F.units() if a % p == 1 and a != 1]
        self.HE = [a for a in F.units() if a % (p ** (n - 1)) == 1 and a != 1]
        # 第 1 跳び i₁（H_K の最小の i(σ) − 1）
        self.i1 = min(F.v(F.sub(F.sigma_pi(a), F.PI)) for a in self.HK) - 1
        self.sig = {a: self._sigma_matrix(a) for a in self.HK}
        self.Ai, self.D = self._basis_inverse()
        self.vD = F.vp(self.D)

    # ---- σ_a の作用行列（π 基底、整数） ----
    def _sigma_matrix(self, a):
        F = self.F
        sp = F.sigma_pi(a)
        cols, term = [], list(F.ONE)
        for _ in range(self.N):
            cols.append(list(term))
            term = F.mul(term, sp)
        return cols

    def _apply(self, cols, x):
        acc = [0] * self.N
        for j in range(self.N):
            xj = x[j]
            if xj == 0:
                continue
            cj = cols[j]
            for r in range(self.N):
                acc[r] += xj * cj[r]
        return acc

    # ---- 基底 {μ^s π^i} → π 基底 の逆行列（整数化） ----
    def _basis_cols(self):
        F = self.F
        cols = []
        for s in range(self.S):
            ms = F.powm(self.mu, s)
            for i in range(self.p):
                cols.append(F.mul(ms, F.powm(F.PI, i)))
        assert len(cols) == self.N
        return cols

    def _basis_inverse(self):
        cols = self._basis_cols()
        n = self.N
        A = [[Fraction(cols[k][r]) for k in range(n)] for r in range(n)]
        I = [[Fraction(1 if i == j else 0) for j in range(n)] for i in range(n)]
        for col in range(n):
            piv = None
            for r in range(col, n):
                if A[r][col] != 0:
                    piv = r
                    break
            assert piv is not None, "singular"
            A[col], A[piv] = A[piv], A[col]
            I[col], I[piv] = I[piv], I[col]
            pv = A[col][col]
            A[col] = [v / pv for v in A[col]]
            I[col] = [v / pv for v in I[col]]
            for r in range(n):
                if r != col and A[r][col] != 0:
                    f = A[r][col]
                    A[r] = [u - f * v for u, v in zip(A[r], A[col])]
                    I[r] = [u - f * v for u, v in zip(I[r], I[col])]
        D = 1
        for row in I:
            for v in row:
                D = D * v.denominator // _gcd(D, v.denominator)
        return [[int(v * D) for v in row] for row in I], D

    def coords(self, y):
        out = []
        for r in range(self.N):
            row = self.Ai[r]
            s = 0
            for j in range(self.N):
                yj = y[j]
                if yj:
                    s += row[j] * yj
            out.append(s)
        return out

    def dist_to_E1(self, y):
        c = self.coords(y)
        out = []
        for i in range(1, self.p):
            best = None
            for s in range(self.S):
                cc = c[s * self.p + i]
                if cc == 0:
                    continue
                w = self.N * (self.F.vp(cc) - self.vD) + self.p * s
                if best is None or w < best:
                    best = w
            if best is not None:
                out.append(best + i)
        return min(out) if out else None

    def eps(self, y):
        best = None
        for a in self.HK:
            v = self.F.v([u - w for u, w in zip(self._apply(self.sig[a], y), y)])
            if v is None:
                continue
            if best is None or v < best:
                best = v
        return best

    def loss(self, y):
        e, d = self.eps(y), self.dist_to_E1(y)
        if e is None or d is None:
            return None
        return e - d


def run(p, n, trials, seed, mod=None):
    t0 = time.time()
    Lr = Layer(p, n)
    mod = mod or p ** 4
    random.seed(seed)
    hist = {}
    worst = None
    for _ in range(trials):
        y = [random.randrange(0, mod) for _ in range(Lr.N)]
        l = Lr.loss(y)
        if l is None:
            continue
        hist[l] = hist.get(l, 0) + 1
        if worst is None or l > worst[0]:
            worst = (l, Lr.eps(y), Lr.dist_to_E1(y), list(y))
    tot = sum(hist.values())
    ok = sum(v for k, v in hist.items() if k <= Lr.i1)
    print(f"=== L = Q_{p}(zeta_{p}^{n})   e = {Lr.N}   i1 = {Lr.i1}   "
          f"|H_K| = {len(Lr.HK)+1} ===")
    print(f"  乱択 {tot} 件（係数 mod {mod}, seed {seed}）")
    print("  損失指数の分布:", dict(sorted(hist.items())))
    print(f"  ★hform（損失 <= i1 = {Lr.i1}）が成り立つ割合 : "
          f"{ok}/{tot} = {100.0*ok/tot:.2f}%")
    print(f"  ★最大損失 = {worst[0]}  (eps={worst[1]}, d={worst[2]})")
    print(f"  ★鋭い一段定数の候補 : c = {p}^({worst[0]}/{Lr.N})"
          f"   /  i1 の上界は {p}^({Lr.i1}/{Lr.N})")
    print(f"  経過 {time.time()-t0:.1f} 秒\n")
    return Lr, worst[0]


def main():
    tri3 = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
    _, m2 = run(3, 3, 20000, 20260909)
    _, m3 = run(3, 4, tri3, 20260909)
    print("=== 法則があるか ===")
    print(f"  k = 2 (e = 18): 最大損失 {m2}  ⇒ 3^({m2}/18)")
    print(f"  k = 3 (e = 54): 最大損失 {m3}  ⇒ 3^({m3}/54)")
    print(f"  axDecay 3 2 = 3^(3/18),  axDecay 3 3 = 3^(3/54)")
    print(f"  ★正規化した指数: k=2 → {Fraction(m2,18)},  k=3 → {Fraction(m3,54)}")
    print(f"  ★axDecay の正規化指数: k=2 → {Fraction(1,6)},  k=3 → {Fraction(1,18)}")


if __name__ == "__main__":
    main()
