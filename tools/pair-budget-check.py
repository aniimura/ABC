# -*- coding: utf-8 -*-
"""厳密整数演算による DeepDescentPairDirect.lean:33-45 の「対の予算」の検算。

★同 :37 は `.../scratchpad/pair/pairsearch.py` を名指ししているが、
★スクリプトはリポジトリに無く、★件数（49,000 / 95,298）は**乱択で seed が記録されていない**
ので★**同じ件数は原理的に再現できない**。★再現できるのは**結論**の方である。

検算対象（同 :39-45 の表、K = ℚ₃(ζ₃)、M = ℚ₃(ζ₂₇)、e_M = 18、深さ k = 2）:

| 予算 | 着地 | v_M での予算 | 測った最悪の必要量 | slack |
|---|---|---|---|---|
| 1 段 axDecay 3 2 | 深さ ≤ 1 | 3 | ★4 | ★−1（届かない） |
| ★対 axDecay 3 1 · axDecay 3 2 | 深さ 0（= K） | ★12 | 8 | ★+4 |

定義（v_M の目盛り）:
  ε      := min_{σ ∈ Gal(M/K), σ≠1} v_M(σx − x)
  1 段の必要量 := ε − v_M(d(x, E₁))     （E₁ = ℚ₃(ζ₉)、深さ ≤ 1 の着地先）
  対の必要量   := ε − v_M(d(x, K))
  1 段の slack := 3  − 1 段の必要量      （axDecay 3 2 = 3^{3/18}）
  対の slack   := 12 − 対の必要量        （axDecay 3 1 · axDecay 3 2 = 3^{12/18}）
"""

import random
import importlib.util
import sys

SRC = r"D:/Math_ABC3/tools/zeta27-distance-check.py"
spec = importlib.util.spec_from_file_location("z27b", SRC)
z = importlib.util.module_from_spec(spec)
sys.modules["z27b"] = z
spec.loader.exec_module(z)

HK = [a for a in range(1, 27) if a % 3 == 1 and a != 1]   # Gal(M/K) \ {1}


def data(x):
    eps = min(z.vM(z.sub(z.apply_sigma(a, x), x)) for a in HK)
    dE = z.dist_to_F(x)   # d(x, Q3(zeta9))
    dK = z.dist_to_K(x)   # d(x, Q3(zeta3))
    if eps is None or dE is None or dK is None:
        return None
    need1 = eps - dE
    need2 = eps - dK
    return eps, dE, dK, need1, need2, 3 - need1, 12 - need2


def main():
    print("=== 前々波で回復した x（WildDescentDistanceOnly の反証候補）===")
    x0 = [71, 78, 9, 50, 30, 54, 67, 43, 40, 3, 25, 11, 40, 44, 74, 63, 38, 24]
    eps, dE, dK, n1, n2, s1, s2 = data(x0)
    print(f"  eps = {eps},  d(x,E1) = {dE},  d(x,K) = {dK}")
    print(f"  1 段: 必要量 = {n1}, 予算 3  ⇒ slack = {s1}   (木の字面 -1)")
    print(f"  対  : 必要量 = {n2}, 予算 12 ⇒ slack = {s2}   (木の字面 +4)")
    print()

    random.seed(20260909)
    trials = 20000
    min1 = None
    min2 = None
    worst1 = None
    hist = {}
    for _ in range(trials):
        y = [random.randrange(0, 81) for _ in range(18)]
        d = data(y)
        if d is None:
            continue
        _, _, _, n1, n2, s1, s2 = d
        min1 = s1 if min1 is None else min(min1, s1)
        min2 = s2 if min2 is None else min(min2, s2)
        worst1 = n1 if worst1 is None else max(worst1, n1)
        hist[(s1, s2)] = hist.get((s1, s2), 0) + 1
    print(f"=== 新しい乱択 {trials} 件（seed = 20260909、★件数は木と違うが結論を測る）===")
    print(f"  ★1 段の最小 slack = {min1}   (木の字面 -1)")
    print(f"  ★対  の最小 slack = {min2}   (木の字面 +4)")
    print(f"  1 段の最悪の必要量 = {worst1}   (木の字面 4)")
    print("  (1 段 slack, 対 slack) の分布（上位 8）:")
    items = sorted(hist.items(), key=lambda kv: -kv[1])
    for k, v in items[:8]:
        print("    ", k, v)
    print()
    bad = [k for k in hist if k[1] < 0]
    print("★対の予算が負になる組み合わせ:", bad if bad else "なし（＝反例ゼロ）")


if __name__ == "__main__":
    main()
