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


def rand_elt(L, tail_only, biased=True):
    """𝒪_L の元を {μ^s π^i} 基底で乱択。tail_only なら i = 0 の成分を 0 にする。"""
    acc = [0] * L.N
    lo = 1 if tail_only else 0
    for s in range(L.S):
        ms = L.F.powm(L.mu, s)
        for i in range(lo, L.p):
            if biased:
                c = L.p ** random.randrange(0, 5) * random.randrange(0, L.p ** 4)
            else:
                c = random.randrange(0, L.p ** 4)
            if c == 0:
                continue
            t = L.F.mul(ms, L.F.powm(L.F.PI, i))
            acc = [a + c * b for a, b in zip(acc, t)]
    return acc


def tail_only_test(p, n, trials, seed=20260909):
    """★★`f₀ = 0`（尾だけ）で測る —— 打ち消しの原因が `f₀` かどうかが 1 回で決まる。

    `f₀ = 0` なら最近点は `x′ = 0` なので `d(x,E₁) = v_L(x)`、
    かつ `‖σx′ − x′‖ = 0 < ‖σy − y‖` は自明に成り立つ
    （`FirstJumpWitness.norm_smul_sub_self_eq_of_lt` の仮説）。
    ⇒ ★もし尾だけで打ち消しが起きなければ、原因は `f₀` だと確定する。
    """
    L = hf.Layer(p, n)
    random.seed(seed)
    iv = {a: L.F.v(L.F.sub(L.F.sigma_pi(a), L.F.PI)) for a in L.HK}
    firstj = [a for a in L.HK if iv[a] == L.i1 + 1]
    for tail_only in (True, False):
        hist = {}
        fj = {}
        dcheck = 0
        for _ in range(trials):
            y = rand_elt(L, tail_only)
            d = L.dist_to_E1(y)
            e = L.eps(y)
            if d is None or e is None:
                continue
            if tail_only and L.F.v(y) == d:
                dcheck += 1
            hist[e - d] = hist.get(e - d, 0) + 1
            best = None
            for a in firstj:
                w = L.F.v([u - v for u, v in zip(sigma_apply(L, a, y), y)])
                if w is None:
                    continue
                if best is None or w < best:
                    best = w
            if best is not None:
                fj[best - d] = fj.get(best - d, 0) + 1
        tot = sum(hist.values())
        tag = "★尾だけ (f0 = 0)" if tail_only else "  対照 (f0 あり)"
        print(f"  {tag}  p={p} n={n} e={L.N} i1={L.i1}  n={tot}")
        print(f"      loss = eps − d の分布       = {dict(sorted(hist.items()))}"
              f"   最大 {max(hist)}")
        print(f"      第1跳びの打ち消し深さの分布 = {dict(sorted(fj.items()))}"
              f"   最大 {max(fj)}")
        if tail_only:
            print(f"      （d = v_L(x) の確認: {dcheck}/{tot}）")


def depth_dense(p, n, want, seed=20260909, maxdraw=400000):
    """★★`d ≡ 1 (mod p)` の点だけを集めて**打ち消しの深さ**の分布を見る。

    `D = d + i₁ = d + (p−1)` を「主項の位」とする。前波で証明した必要条件は
    `p ∣ D`（⟺ `d ≡ 1 mod p`）だった。★母集団をそこに絞ると打ち消しが `p` 倍濃く出る。

    深さ = `min_{σ: i(σ)=i₁+1} v_L(σx − x) − D`。
    ★同時に、その最小を実現する **𝒪_{E₁} 成分の番号 j** も記録する
    （位は `j (mod p)` に合同なので、`D` の上の「空き枠」は
      `D+1, …, D+(p−1)`（尾の成分）と `D+p`（`E₁` 成分）である）。
    """
    L = hf.Layer(p, n)
    random.seed(seed)
    iv = {a: L.F.v(L.F.sub(L.F.sigma_pi(a), L.F.PI)) for a in L.HK}
    firstj = [a for a in L.HK if iv[a] == L.i1 + 1]
    depth = {}
    slot = {}
    got = 0
    draws = 0
    while got < want and draws < maxdraw:
        draws += 1
        y = rand_elt(L, False)
        d = L.dist_to_E1(y)
        if d is None or (d - 1) % p != 0:
            continue
        got += 1
        best = None
        besta = None
        for a in firstj:
            w = L.F.v([u - v for u, v in zip(sigma_apply(L, a, y), y)])
            if w is None:
                continue
            if best is None or w < best:
                best, besta = w, a
        if best is None:
            continue
        D = d + L.i1
        depth[best - D] = depth.get(best - D, 0) + 1
        Av = coords_v(L, [u - v for u, v in
                          zip(sigma_apply(L, besta, y), y)])
        j = min((jj for jj in range(p) if Av[jj] is not None),
                key=lambda jj: Av[jj] + jj)
        slot[(best - D, j)] = slot.get((best - D, j), 0) + 1
    print(f"  p={p} n={n} e={L.N} i1={L.i1}   d≡1 (mod {p}) の点 {got} 件"
          f"（抽選 {draws} 回）")
    print(f"    ★深さ = v(σx−x) − D の分布 = {dict(sorted(depth.items()))}"
          f"   ★最大 {max(depth) if depth else None}")
    print(f"    深さと成分 j の対応 = "
          f"{dict(sorted((k, v) for k, v in slot.items()))}")


def dense():
    print("=== ★★d ≡ 1 (mod p) に絞った打ち消しの深さ ===")
    depth_dense(3, 3, 4000)
    depth_dense(2, 4, 4000)
    depth_dense(3, 4, 1200)
    depth_dense(2, 5, 2000)
    depth_dense(5, 3, 250)


def tail():
    print("=== ★★尾だけ（f₀ = 0）で打ち消しは起きるか ===")
    tail_only_test(3, 3, 6000)
    tail_only_test(3, 4, 1500)
    tail_only_test(2, 4, 6000)
    tail_only_test(2, 5, 3000)


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
