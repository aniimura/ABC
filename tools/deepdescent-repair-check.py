# -*- coding: utf-8 -*-
"""厳密整数演算による DeepDescentRepair.lean:29-45 の反例の再導出。

★同ファイル :110 は「再現スクリプトは …/scratchpad/…」と書いており、
★スクリプトはリポジトリに無い。★ただし x は明示されている（:35）ので再導出できる。

設定（DeepDescentRepair.lean:31-35 の字面）:
  F = ℚ₃(ζ₃)、L = ℚ₃(ζ₂₇)（e_L = 18、Gal(L/F) = ℤ/9 巡回）、
  π = ζ₂₇ − 1、π_E = ζ₉ − 1（v_L(π_E) = 3）、E₁ = ℚ₃(ζ₉)。
  ★x := π³ + 2·π_E + π_E⁵

検算対象（同 :37-42 の表）:
  [F(x):F] = 9 （⇒ wildDepth = 2）
  min_{σ≠1} v_L(σx − x) = 23   （8 元すべて: 23,23,27,23,23,27,23,23）
  d(x, E₁) = 19
  d(x, F)  = 15
"""

import importlib.util
import sys

SRC = r"D:/Math_ABC3/tools/zeta27-distance-check.py"
spec = importlib.util.spec_from_file_location("z27", SRC)
z = importlib.util.module_from_spec(spec)
sys.modules["z27"] = z
spec.loader.exec_module(z)


def main():
    PI = z.PI
    ONE = z.ONE
    mul, sub, vM, powm = z.mul, z.sub, z.vM, z.pow_mod
    add = lambda a, b: [u + v for u, v in zip(a, b)]

    # π_E = ζ₉ − 1 = (1+π)^3 − 1   （z の MU がこれ）
    PI_E = z.MU
    print("v_L(pi)   =", vM(PI))
    print("v_L(pi_E) =", vM(PI_E), "  (期待 3)")

    # x = π³ + 2·π_E + π_E⁵
    x = add(add(powm(PI, 3), [2 * c for c in PI_E]), powm(PI_E, 5))
    print("v_L(x)    =", vM(x))
    print()

    # Gal(L/F), F = Q3(zeta3):  a = 1 mod 3
    HF = [a for a in range(1, 27) if a % 3 == 1]
    vals = []
    orbit = set()
    for a in HF:
        y = z.apply_sigma(a, x)
        orbit.add(tuple(y))
        if a != 1:
            vals.append((a, vM(sub(y, x))))
    print("Gal(L/F) の元 a と v_L(sigma_a x - x):")
    for a, v in vals:
        print(f"   a={a:3d}  v = {v}")
    print("   ★min =", min(v for _, v in vals), "  (期待 23)")
    print("   ★軌道の大きさ [F(x):F] =", len(orbit), "  (期待 9)")
    print()
    print("d(x, E1 = Q3(zeta9)) =", z.dist_to_F(x), "  (期待 19)")
    print("d(x, F  = Q3(zeta3)) =", z.dist_to_K(x), "  (期待 15)")
    print()
    eps = min(v for _, v in vals)
    need = eps - 3          # axDecay 3 2 = 3^{3/18} ⇒ 要求は v = eps - 3
    print(f"配られた字面の要求: v_L(x - x') >= {need}   (eps = {eps}, axDecay 3 2 = 3^(3/18))")
    dE, dF = z.dist_to_F(x), z.dist_to_K(x)
    print(f"  E1 へ: v = {dE}  >= {need} ?  ->  {dE >= need}   (★満たないなら証人無し)")
    print(f"  F  へ: v = {dF}  >= {need} ?  ->  {dF >= need}")
    print()
    if dE < need and dF < need:
        print("★★ L の中に証人は無い —— DeepDescentRepair.lean:44 と一致")


if __name__ == "__main__":
    main()
