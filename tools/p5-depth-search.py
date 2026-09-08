# -*- coding: utf-8 -*-
"""★`p = 5` で打ち消しの深さ 2 を出す（＝ `loss ≤ p+1` の真偽を決める）。

★盲目的な乱択だと 10^4〜10^5 件（10〜100 分）かかると見積もった。
★本スクリプトは**層別に構成する**ことで 10^3 件（数分）に落とす。

構成:
  x = Σ_{i<p} μ^{s_i} u_i π^i   （μ は E₁ の素元、v_L(μ) = p、u_i は単数）
  ⇒ v_L(f_i) = p·s_i,  d = min_{i≥1}(p·s_i + i),  D = d + (p−1)

打ち消し（深さ ≥ 1）が起きる条件は  v_L(σ f₀ − f₀) = D。
`t := v_L(σμ − μ) − v_L(μ)` を測れば  v_L(σ f₀ − f₀) = p·s₀ + t  なので
  ★p·s₀ + t = D  から s₀ が決まる。
⇒ ★s₀ を当てて u₀ を振れば、打ち消しが**設計で**起きる。

深さ 2 はさらに「B₁ が最小の枠を外す」ことを要する（前波の SlotStructure）。
★これは u_i を振って探す。
"""

import importlib.util
import random
import sys
import time

K3 = r"D:/Math_ABC3/tools/hform-k3-check.py"
_spec = importlib.util.spec_from_file_location("hfk3c", K3)
hf = importlib.util.module_from_spec(_spec)
sys.modules["hfk3c"] = hf
_spec.loader.exec_module(hf)

C = r"D:/Math_ABC3/tools/cancel-structure-check.py"
_spec2 = importlib.util.spec_from_file_location("csc", C)
cs = importlib.util.module_from_spec(_spec2)
sys.modules["csc"] = cs
_spec2.loader.exec_module(cs)


def build(L, expo, units):
    """x = Σ_i μ^{expo[i]} · units[i] · π^i。expo[i] が None ならその項なし。"""
    acc = [0] * L.N
    for i in range(L.p):
        if expo[i] is None:
            continue
        term = L.F.mul(L.F.powm(L.mu, expo[i]), L.F.powm(L.F.PI, i))
        c = units[i]
        acc = [a + c * b for a, b in zip(acc, term)]
    return acc


def measure_t(L, a):
    """t = v_L(σ_a μ − μ) − v_L(μ)。"""
    sp = L.sig[a]
    smu = L._apply(sp, L.mu)
    return L.F.v(L.F.sub(smu, L.mu)) - L.F.v(L.mu)


def depth_of(L, y, firstj):
    d = L.dist_to_E1(y)
    if d is None:
        return None, None, None
    best, besta = None, None
    for a in firstj:
        w = L.F.v([u - v for u, v in zip(L._apply(L.sig[a], y), y)])
        if w is None:
            continue
        if best is None or w < best:
            best, besta = w, a
    if best is None:
        return None, None, None
    D = d + L.i1
    return best - D, d, besta


def run(p, n, trials, s1=2, seed=20260909, verbose=True):
    t0 = time.time()
    L = hf.Layer(p, n)
    iv = {a: L.F.v(L.F.sub(L.F.sigma_pi(a), L.F.PI)) for a in L.HK}
    firstj = [a for a in L.HK if iv[a] == L.i1 + 1]
    t = measure_t(L, firstj[0])
    if verbose:
        print(f"=== p={p} n={n} e={L.N} i1={L.i1}  第1跳びの元 {len(firstj)} 個 ===")
        print(f"  v_L(mu) = {L.F.v(L.mu)},  t = v_L(sigma mu - mu) - v_L(mu) = {t}")
    # d は i = 1 で達成させる: s_i (i>=2) を大きめに取る
    d = p * s1 + 1
    D = d + L.i1
    if (D - t) % p != 0:
        print(f"  ★D - t = {D - t} が p で割り切れない ⇒ この s1 では打ち消し不能")
        return None
    s0 = (D - t) // p
    if s0 < 0:
        print(f"  ★s0 = {s0} < 0 ⇒ この s1 では不可。s1 を上げる。")
        return None
    if verbose:
        print(f"  s1 = {s1} ⇒ d = {d}, D = {D}, ★s0 = {s0}"
              f"（設計で v(sigma f0 - f0) = D）")
    random.seed(seed)
    hist = {}
    found = None
    for _ in range(trials):
        # ★観測された最悪点は v(f_1) = v(f_2) = … だった（SlotStructure の表）
        expo = [s0] + [s1] * (p - 1)
        units = [random.randrange(1, p ** 4) for _ in range(p)]
        for i in range(p):
            if units[i] % p == 0:
                units[i] += 1
        y = build(L, expo, units)
        dep, dd, _ = depth_of(L, y, firstj)
        if dep is None or dd != d:
            continue
        hist[dep] = hist.get(dep, 0) + 1
        if dep >= 2 and found is None:
            found = (dep, list(units))
    if verbose:
        print(f"  ★深さの分布 = {dict(sorted(hist.items()))}"
              f"   最大 {max(hist) if hist else None}   ({time.time()-t0:.1f}s)")
        if found:
            print(f"  ★★★深さ {found[0]} を発見。units = {found[1]}")
    return hist


def main():
    tri = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    print("### 対照: p = 3（既知の答え 最大 2）###")
    for s1 in (1, 2, 3):
        run(3, 3, tri, s1=s1)
    print()
    print("### ★本番: p = 5 ###")
    for s1 in (1, 2, 3):
        run(5, 3, tri, s1=s1)


if __name__ == "__main__":
    main()
