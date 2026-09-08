# -*- coding: utf-8 -*-
"""★「打ち消しはなぜ 2 目盛りで止まるのか」を、最悪点の構造から読む。

前波までに測ったこと（tools/hform-general-p-check.py）:
  損失 = eps − d(x,E₁) は 層の付値の目盛りで測ると p に依らず
    生成的 = i₁ = p−1、最大 = p+1（p = 2,3 で確認）。
★なぜ最大が p+1 で止まるのかは分かっていない。

本スクリプトは L = ℚ₃(ζ₂₇)（e = 18）で損失 4 の x を集め、
𝒪_{E₁} 係数 x = f₀ + f₁π + f₂π² と σx − x = A₀ + A₁π + A₂π² を厳密に展開して、
★どの成分がどの位で打ち消されているかを見る。

読み方:
  v_L(f_i) は 3ℤ（E₁ の付値群）。
  d = min_i (v(f_i) + i)  (i = 1,2)
  σ ∉ H_E なら v(σπ − π) = 3（第 1 跳び）、σ ∈ H_E なら 9。
  生成的には v(σx − x) = v(f₁) + 3 = d + 2。
"""

import importlib.util
import random
import sys

K3 = r"D:/Math_ABC3/tools/hform-k3-check.py"
_spec = importlib.util.spec_from_file_location("hfk3b", K3)
hf = importlib.util.module_from_spec(_spec)
sys.modules["hfk3b"] = hf
_spec.loader.exec_module(hf)


def vE(L, c):
    """整数化された座標 c（D で割る前）の v_L。None なら 0。"""
    if c == 0:
        return None
    return L.N * (L.F.vp(c) - L.vD)


def coords_v(L, y):
    """y の 𝒪_{E₁} 係数 f_0..f_{p-1} の v_L を返す（None は 0）。"""
    c = L.coords(y)
    out = []
    for i in range(L.p):
        best = None
        for s in range(L.S):
            cc = c[s * L.p + i]
            if cc == 0:
                continue
            w = L.N * (L.F.vp(cc) - L.vD) + L.p * s
            if best is None or w < best:
                best = w
        out.append(best)
    return out


def sigma_apply(L, a, y):
    return L._apply(L.sig[a], y)


def report_point(L, y, tag):
    p = L.p
    fv = coords_v(L, y)
    d = L.dist_to_E1(y)
    e = L.eps(y)
    print(f"  [{tag}] v(f_i) = {fv}   d = {d}   eps = {e}   loss = {e - d}")
    rows = []
    for a in L.HK:
        iv = L.F.v(L.F.sub(L.F.sigma_pi(a), L.F.PI))
        w = L.F.v([u - v for u, v in zip(sigma_apply(L, a, y), y)])
        rows.append((a, iv, w))
    rows.sort(key=lambda r: r[2] if r[2] is not None else 10 ** 9)
    for a, iv, w in rows[:4]:
        Av = coords_v(L, [u - v for u, v in zip(sigma_apply(L, a, y), y)])
        pred = None if fv[1] is None else fv[1] + iv
        print(f"      sigma_{a:2d} i={iv}  v(sx-x)={w}"
              f"   A 成分 v(A_j)+j = "
              f"{[None if Av[j] is None else Av[j]+j for j in range(p)]}"
              f"   生成的な予測 v(f1)+i = {pred}"
              f"   ずれ {None if (w is None or pred is None) else w - pred}")


def main():
    L = hf.Layer(3, 3)
    print(f"=== L = Q3(zeta27)  e = {L.N}  i1 = {L.i1} ===")
    print("  H_E（E₁ を固定）=", L.HE)
    print()
    random.seed(20260909)
    worst, generic, mid = [], [], []
    for _ in range(60000):
        y = [random.randrange(0, 81) for _ in range(L.N)]
        l = L.loss(y)
        if l is None:
            continue
        if l == 4 and len(worst) < 6:
            worst.append(list(y))
        elif l == 3 and len(mid) < 3:
            mid.append(list(y))
        elif l == 2 and len(generic) < 3:
            generic.append(list(y))
        if len(worst) >= 6 and len(mid) >= 3 and len(generic) >= 3:
            break
    print("=== 生成的（損失 2 = i₁）===")
    for i, y in enumerate(generic):
        report_point(L, y, f"gen{i}")
    print()
    print("=== 損失 3（= axDecay）===")
    for i, y in enumerate(mid):
        report_point(L, y, f"mid{i}")
    print()
    print("=== ★損失 4（最大）===")
    for i, y in enumerate(worst):
        report_point(L, y, f"max{i}")


def top_component_test(p, n, trials, seed=20260909, mod=None):
    """★読めた構造の検証 —— 「最上位成分 A_{p−1} だけで上界が出る」か。

    最悪点の観察: 損失 4 の点では常に  v(A_2) + 2 = d + 4  が最小成分だった。
    ⇒ 主張: ある σ ∈ H_K \\ H_E で  v(A_{p−1}) + (p−1) ≤ d + p + 1。
    ★これが正しければ、eps ≤ v(A_{p−1}) + (p−1) （min の性質）から loss ≤ p+1 が出る。
    """
    L = hf.Layer(p, n)
    mod = mod or p ** 4
    random.seed(seed)
    HKonly = [a for a in L.HK if a not in L.HE]
    hist = {}
    bad = 0
    for _ in range(trials):
        y = [p ** random.randrange(0, 5) * random.randrange(0, p ** 4)
             for _ in range(L.N)]
        d = L.dist_to_E1(y)
        if d is None:
            continue
        best = None
        for a in HKonly:
            Av = coords_v(L, [u - v for u, v in
                              zip(sigma_apply(L, a, y), y)])
            t = Av[p - 1]
            if t is None:
                continue
            w = t + (p - 1)
            if best is None or w < best:
                best = w
        if best is None:
            continue
        g = best - d
        hist[g] = hist.get(g, 0) + 1
        if g > p + 1:
            bad += 1
    print(f"  p={p} n={n} e={L.N}: min_σ (v(A_{{p-1}})+(p-1)) − d の分布 = "
          f"{dict(sorted(hist.items()))}")
    print(f"    ★ p+1 = {p+1} を超えた件数 = {bad} / {sum(hist.values())}"
          f"   最大 = {max(hist)}")


def firstjump_cancel_test(p, n, trials, seed=20260909):
    """★★読めた法則の検証 —— **第 1 跳びの元での打ち消しは高々 2 目盛り**。

    σ を i(σ) = i₁ + 1（第 1 跳びの層）の元とすると、生成的には
        v(σx − x) = v(f₁) + i(σ) = d + i₁   （d = v(f₁) + 1 のとき）
    である。★測定した最悪点では、この値が **ちょうど 2 目盛り**だけ深くなっていた。
    ⇒ 主張: min_{σ: i(σ)=i₁+1} v(σx − x) ≤ d + i₁ + 2。
    ★これが正しければ loss ≤ i₁ + 2 = p + 1 が出る（eps はさらに小さいから）。
    """
    L = hf.Layer(p, n)
    random.seed(seed)
    iv = {a: L.F.v(L.F.sub(L.F.sigma_pi(a), L.F.PI)) for a in L.HK}
    first = [a for a in L.HK if iv[a] == L.i1 + 1]
    hist = {}
    bad = 0
    for _ in range(trials):
        y = [p ** random.randrange(0, 5) * random.randrange(0, p ** 4)
             for _ in range(L.N)]
        d = L.dist_to_E1(y)
        if d is None:
            continue
        best = None
        for a in first:
            w = L.F.v([u - v for u, v in zip(sigma_apply(L, a, y), y)])
            if w is None:
                continue
            if best is None or w < best:
                best = w
        if best is None:
            continue
        g = best - d
        hist[g] = hist.get(g, 0) + 1
        if g > L.i1 + 2:
            bad += 1
    tot = sum(hist.values())
    print(f"  p={p} n={n} e={L.N} i1={L.i1}  |第1跳びの元| = {len(first)}")
    print(f"    min_σ v(σx−x) − d の分布 = {dict(sorted(hist.items()))}")
    print(f"    ★ i1+2 = {L.i1+2} を超えた件数 = {bad} / {tot}"
          f"   最大 = {max(hist)}")


def first():
    print("=== ★★第 1 跳びの元での打ち消しは高々 2 目盛りか ===")
    firstjump_cancel_test(3, 3, 8000)
    firstjump_cancel_test(3, 4, 2000)
    firstjump_cancel_test(2, 4, 8000)
    firstjump_cancel_test(2, 5, 4000)
    firstjump_cancel_test(5, 3, 400)


def top():
    print("=== ★最上位成分だけで上界が出るか ===")
    top_component_test(3, 3, 8000)
    top_component_test(3, 4, 2000)
    top_component_test(2, 4, 8000)
    top_component_test(2, 5, 4000)
    top_component_test(5, 2, 3000)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "top":
        top()
    else:
        main()
