# -*- coding: utf-8 -*-
"""厳密整数演算による「正規化した跡の欠損 γ」の検算。

★EquivariantProjectionDescent.lean:19 は
  スクリプト `.../scratchpad/grad/three.py`, `trnorm.py`, `maxloss.py`
を名指ししているが、★リポジトリに無い。本スクリプトはその中心の主張を再導出する。

検算対象（同 :26-31）:
  P = (1/p)·Tr_{M/E} の作用素ノルム `Q = ‖P‖ = p^{γ/e_M}`、
  ★γ = p − 1 が測った 4 本の層すべてで一致した:
    ℚ₃(ζ₂₇)/ℚ₃(ζ₉): 2、ℚ₃(ζ₈₁)/ℚ₃(ζ₂₇): 2、
    ℚ₂(ζ₈)/ℚ₂(ζ₄): 1、ℚ₂(ζ₁₆)/ℚ₂(ζ₈): 1
  および同 :23-25 の t = 8（上の層の跳び）、s = 6、m = 2。

方法:
  M = ℚ_p(ζ_{p^n})、π = ζ_{p^n} − 1、𝒪_M = ℤ_p[π]、v_M(Σ b_i π^i) = min_i(e·v_p(b_i)+i)。
  E = ℚ_p(ζ_{p^{n-1}})、Gal(M/E) = { σ_a : a ≡ 1 mod p^{n-1} }（位数 p）。
  Tr は 𝒪_E 線形で 𝒪_M = 𝒪_E[π]（基底 {1, π, …, π^{p−1}}）。
  ★作用素ノルムは「比の sup」なので
    γ = max_{i<p} ( v_M(π^i) − v_M(P π^i) ) = max_{i<p} ( i + e_M − v_M(Tr(π^i)) )。
  ★最初の実装は min を取って γ = 0 を出した（誤り）。sup と min を取り違えていた。

  ★古典的な裏付け（Serre, Corps Locaux V §3 Lemme 4）:
    Tr(𝒪_M) = 𝔭_E^{floor(d/e)}、d = (p−1)(i+1)、e = e(M/E) = p。
    M/E = ℚ₃(ζ₂₇)/ℚ₃(ζ₉) なら d = 2·9 = 18、floor(18/3) = 6
    ⇒ v_M(Tr 𝒪_M) = 3·6 = 18 = e_M ⇒ γ = max_i i = p−1。
    ★本スクリプトの数値と一致する（2 つの独立な道が合う）。
"""

import importlib.util
import sys

SRC = r"D:/Math_ABC3/tools/zeta-tower-check.py"
spec = importlib.util.spec_from_file_location("ztower", SRC)
zt = importlib.util.module_from_spec(spec)
sys.modules["ztower"] = zt
spec.loader.exec_module(zt)

Field = zt.Field


def trace_defect(p, n):
    """M = Q_p(zeta_{p^n}) / E = Q_p(zeta_{p^{n-1}}) の gamma を返す。"""
    F = Field(p, n)
    sub = p ** (n - 1)
    G = [a for a in F.units() if a % sub == 1]
    assert len(G) == p, (len(G), p)
    vals = []
    for i in range(p):
        y = F.powm(F.PI, i)
        tr = [0] * F.N
        for a in G:
            sy = F.red(subst(F, a, y))
            tr = [u + v for u, v in zip(tr, sy)]
        vals.append((i, F.v(tr)))
    # ★作用素ノルムは「比の sup」であって「v の min」ではない。
    # y = Σ c_i π^i (c_i ∈ ᵒ_E) について ‖y‖ = max_i ‖c_i π^i‖ なので
    #   γ = max_i ( v_M(π^i) - v_M(P π^i) ) = max_i ( i + e_M - v_M(Tr π^i) )
    gamma = max(i + F.N - v for i, v in vals if v is not None)
    return F, G, vals, gamma


def subst(F, a, y):
    """σ_a(y) —— y = Σ b_j π^j に σ_a(π) を代入。"""
    sp = F.sigma_pi(a)
    acc = [0] * F.N
    term = list(F.ONE)
    for j in range(F.N):
        if y[j] != 0:
            acc = [u + y[j] * v for u, v in zip(acc, term)]
        term = F.mul(term, sp)
    return acc


def main():
    cases = [(3, 3), (3, 4), (2, 3), (2, 4)]
    print("M = Q_p(zeta_{p^n})  ⊃  E = Q_p(zeta_{p^{n-1}})")
    print()
    for p, n in cases:
        F, G, vals, gamma = trace_defect(p, n)
        print(f"p={p}, n={n}:  e_M = {F.N},  Gal(M/E) = {sorted(G)}")
        print("   v_M(Tr(pi^i)) :", {i: v for i, v in vals})
        print(f"   ★gamma = max_i (i + e_M - v_M(Tr pi^i)) = {gamma}"
              f"   (期待 p-1 = {p-1})   ->  {gamma == p - 1}")
        print()

    # 測定 1 のその他の数値（p=3, n=3）
    F = Field(3, 3)
    ivals = {a: F.v(F.sub(F.sigma_pi(a), F.PI)) for a in F.units()}
    subE = [a for a in F.units() if a % 9 == 1 and a != 1]
    subK = [a for a in F.units() if a % 3 == 1 and a != 1]
    t = min(ivals[a] for a in subE) - 1
    m = min(ivals[a] for a in subK) - 1
    print("=== EquivariantProjectionDescent.lean:23-25 の t / s / m ===")
    print(f"  t = 上の層 M/E = Q3(zeta27)/Q3(zeta9) の跳び = {t}   (期待 8)")
    print(f"  m = 第 1 跳び(= Gal(M/K) の最小の伸び − 1) = {m}   (期待 2)")
    # s: 下の層 E/K の跳び 2 を v_M に換算 = 2 * e(M/E) = 2*3 = 6
    FE = Field(3, 2)
    ivE = {a: FE.v(FE.sub(FE.sigma_pi(a), FE.PI)) for a in FE.units()}
    subEK = [a for a in FE.units() if a % 3 == 1 and a != 1]
    sE = (min(ivE[a] for a in subEK) - 1) if subEK else None
    print(f"  下の層 E/K = Q3(zeta9)/Q3(zeta3) の跳び(E の目盛り) = {sE}   (期待 2)")
    print(f"  s = その v_M 換算 = {sE} * e(M/E) = {sE} * 3 = {sE * 3 if sE is not None else None}"
          f"   (期待 6)")
    if sE is not None:
        gamma = 2
        print()
        print("=== 測定 2 の合成（p=3, n=3, 底 K = Q3(zeta3)）===")
        print(f"  素朴な telescope t + s          = {t + sE*3}   (期待 14、予算 12 を超える)")
        print(f"  収縮つき max(t, s, t+s-m)       = "
              f"{max(t, sE*3, t + sE*3 - m)}   (期待 12)")
        print(f"  ★射影つき max(t+gamma, s+gamma) = "
              f"{max(t + gamma, sE*3 + gamma)}   (期待 10)")


if __name__ == "__main__":
    main()
