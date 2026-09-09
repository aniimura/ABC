# -*- coding: utf-8 -*-
"""★★手で導いた成分の公式を**機械で確かめる**（本日まだ一度も検算していない）。

`TopIndexSurvives.lean` の 5 段の (2) で私はこう書いた:

  `f_i·((π+ρ)^i − π^i)` の主項は `f_i·i·π^{i−1}·ρ = i f_i w (π^{i−1} + π^i)`
  ⇒ ★`B_j = w·( j·f_j + (j+1)·f_{j+1} ) + (より深い項)`

★これは**手計算であって、数値で確かめていない**。本スクリプトが確かめる。

やること: `x` を乱択し、`𝒪_{E₁}` 係数 `f_i` を厳密に取り出し、
予測 `P_j := w·( j·f_j + (j+1)·f_{j+1} )` と実際の `B_j`（`σx − x` の成分）を比べ、
★`v(B_j − P_j) > v(P_j)`（＝ 予測が主部）を検査する。
"""

import importlib.util
import random
import sys

K3 = r"D:/Math_ABC3/tools/hform-k3-check.py"
_spec = importlib.util.spec_from_file_location("hfk3d", K3)
hf = importlib.util.module_from_spec(_spec)
sys.modules["hfk3d"] = hf
_spec.loader.exec_module(hf)


def elt_from_coords(L, c, i):
    """成分 i の係数（D 倍された `𝒪_{E₁}` の元）を L の元として組み立てる。"""
    acc = [0] * L.N
    for s in range(L.S):
        v = c[s * L.p + i]
        if v == 0:
            continue
        ms = L.F.powm(L.mu, s)
        acc = [a + v * b for a, b in zip(acc, ms)]
    return acc


def add(L, a, b):
    return [x + y for x, y in zip(a, b)]


def smul(L, k, a):
    return [k * x for x in a]


def run(p, n, trials, seed=20260909):
    L = hf.Layer(p, n)
    F = L.F
    iv = {a: F.v(F.sub(F.sigma_pi(a), F.PI)) for a in L.HK}
    firstj = [a for a in L.HK if iv[a] == L.i1 + 1]
    random.seed(seed)
    ok = bad = skip = 0
    examples = []
    for _ in range(trials):
        y = [p ** random.randrange(0, 4) * random.randrange(0, p ** 3)
             for _ in range(L.N)]
        a = firstj[0]
        b = (a - 1) // p
        # w = zeta_{p^{n-1}}^b - 1
        one_plus_pi = F.red([1, 1] + [0] * (F.N - 2))
        w = F.sub(F.powm(one_plus_pi, p * b), F.ONE)
        cx = L.coords(y)
        Df = [elt_from_coords(L, cx, i) for i in range(p)]
        Df.append([0] * L.N)                       # f_p = 0
        diff = [u - v for u, v in zip(L._apply(L.sig[a], y), y)]
        cb = L.coords(diff)
        for j in range(p):
            DB = elt_from_coords(L, cb, j)
            inner = add(L, smul(L, j, Df[j]), smul(L, j + 1, Df[j + 1]))
            DP = F.mul(w, inner)
            vP = F.v(DP)
            if vP is None:
                skip += 1
                continue
            vR = F.v(F.sub(DB, DP))
            if vR is None or vR > vP:
                ok += 1
            else:
                bad += 1
                if len(examples) < 3:
                    examples.append((j, vP, vR))
    print(f"  p={p} n={n}: 成分の検査 {ok+bad} 件（P=0 で飛ばした {skip} 件）")
    print(f"    ★v(B_j − P_j) > v(P_j) が成立 {ok} 件 / ★不成立 {bad} 件")
    if examples:
        print(f"    反例（j, v(P), v(B−P)）: {examples}")


def main():
    print("=== 手で導いた B_j = w(j f_j + (j+1) f_{j+1}) + 深い項 の検算 ===")
    run(3, 3, 200)
    run(3, 4, 60)
    run(5, 3, 20)
    run(2, 4, 200)


if __name__ == "__main__":
    main()
