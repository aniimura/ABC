# -*- coding: utf-8 -*-
"""★`ComponentExpansion.lean` が言っていることを塔で測る（規律 9: 式には母集団）。

Lean 側は「凍結模型」を証明した:

  σx − x = Σ_j f_j·((π+ρ)^j − π^j)      ← ★σ が f_j を動かさないと仮定した形
         = w·Σ_j ( j f_j + (j+1) f_{j+1} )·π^j + E,   ‖E‖ ≤ ‖ρ‖²
  （ρ = w(1+π), f_p = 0, ‖f_j‖ ≤ 1, ‖π‖ ≤ 1, ‖ρ‖ ≤ 1）

塔では σ ∈ H_K（a ≡ 1 mod p）は E₁ を**動かす**（a ≡ 1 mod p^{n-1} でない限り）。
そこで測るのは 5 点:

  (m1) ρ = w·(1+π) が厳密に成り立つか（Lean の仮定 hw）
  (m2) v( 凍結模型 − 主部 ) ≥ 2·v(ρ)          ← ★Lean の定理そのもの
  (m3) v( 実際の σx−x − 凍結模型 ) ≥ 2·v(ρ)   ← ★模型と塔のずれ（σ が f_j を動かす分）
  (m4) v( 実際の σx−x − 主部 ) ≥ 2·v(ρ)       ← ★結論（塔でも成り立つか）
  (m5) σ(f_j) = f_j か（★動かすなら模型は塔と一致しない。何件動くかを数える）

すべて D 倍した `𝒪_{E₁}` 係数で計算する（D は有理整数なので σ と可換、付値の比較に影響しない）。
"""

import importlib.util
import random
import sys

K3 = r"D:/Math_ABC3/tools/hform-k3-check.py"
_spec = importlib.util.spec_from_file_location("hfk3exp", K3)
hf = importlib.util.module_from_spec(_spec)
sys.modules["hfk3exp"] = hf
_spec.loader.exec_module(hf)


def elt_from_coords(L, c, i):
    acc = [0] * L.N
    for s in range(L.S):
        v = c[s * L.p + i]
        if v == 0:
            continue
        ms = L.F.powm(L.mu, s)
        acc = [a + v * b for a, b in zip(acc, ms)]
    return acc


def add(a, b):
    return [x + y for x, y in zip(a, b)]


def sub(a, b):
    return [x - y for x, y in zip(a, b)]


def smul(k, a):
    return [k * x for x in a]


def vge(F, z, bound):
    """v(z) ≥ bound か（z = 0 は ∞ とみなす）。"""
    v = F.v(z)
    return (v is None) or (v >= bound), v


def run(p, n, trials, seed=20260909):
    L = hf.Layer(p, n)
    F = L.F
    iv = {a: F.v(F.sub(F.sigma_pi(a), F.PI)) for a in L.HK}
    i1 = L.i1
    a0 = [a for a in L.HK if iv[a] == i1 + 1][0]
    b = (a0 - 1) // p
    one_plus_pi = F.red([1, 1] + [0] * (F.N - 2))
    w = F.sub(F.powm(one_plus_pi, p * b), F.ONE)
    sp = F.sigma_pi(a0)                     # σπ = π + ρ
    rho = F.sub(sp, F.PI)
    vrho = F.v(rho)
    m1 = (F.mul(w, one_plus_pi) == rho)     # (m1)

    random.seed(seed)
    c2 = c3 = c4 = 0
    fixed = moved = 0
    m5b = [0, 0, None]                      # [成立数, 総数, 最小の v]
    worst = {"m2": None, "m3": None, "m4": None}
    for _ in range(trials):
        y = [p ** random.randrange(0, 4) * random.randrange(0, p ** 3)
             for _ in range(L.N)]
        cx = L.coords(y)
        Df = [elt_from_coords(L, cx, i) for i in range(p)]
        Df.append([0] * L.N)                # f_p = 0

        # (m5) σ が成分を動かすか / (m5b) 動く量は 2v(ρ) 以深か
        for j in range(p):
            gj = L._apply(L.sig[a0], Df[j])
            if gj == Df[j]:
                fixed += 1
            else:
                moved += 1
            ok5b, v5b = vge(F, sub(gj, Df[j]), 2 * vrho)
            m5b[0] += ok5b
            m5b[1] += 1
            if v5b is not None and (m5b[2] is None or v5b < m5b[2]):
                m5b[2] = v5b

        # 凍結模型 Σ_j f_j ((π+ρ)^j − π^j)
        frozen = [0] * L.N
        for j in range(p):
            t = F.sub(F.powm(sp, j), F.powm(F.PI, j))
            frozen = add(frozen, F.mul(Df[j], t))

        # 主部 w Σ_j ( j f_j + (j+1) f_{j+1} ) π^j
        main = [0] * L.N
        for j in range(p):
            coef = add(smul(j, Df[j]), smul(j + 1, Df[j + 1]))
            main = add(main, F.mul(coef, F.powm(F.PI, j)))
        main = F.mul(w, main)

        # 実際の D(σx − x)（D は有理整数なので σ と可換）
        actual = [0] * L.N
        for j in range(p):
            actual = add(actual, F.mul(L._apply(L.sig[a0], Df[j]), F.powm(sp, j)))
            actual = sub(actual, F.mul(Df[j], F.powm(F.PI, j)))

        ok2, v2 = vge(F, sub(frozen, main), 2 * vrho)
        ok3, v3 = vge(F, sub(actual, frozen), 2 * vrho)
        ok4, v4 = vge(F, sub(actual, main), 2 * vrho)
        c2 += ok2
        c3 += ok3
        c4 += ok4
        for k, v in (("m2", v2), ("m3", v3), ("m4", v4)):
            if v is not None and (worst[k] is None or v < worst[k]):
                worst[k] = v

    print(f"  p={p} n={n} 標本 {trials}  v(ρ)={vrho}  2v(ρ)={2 * vrho}"
          f"  v(w)={F.v(w)}  i1={i1}")
    print(f"    (m1) ρ = w(1+π)              : {'成立' if m1 else '★不成立'}")
    print(f"    (m2) 凍結模型 − 主部 が深い  : {c2}/{trials}"
          f"  最小の v = {worst['m2']}")
    print(f"    (m3) 実際 − 凍結模型 が深い  : {c3}/{trials}"
          f"  最小の v = {worst['m3']}")
    print(f"    (m4) 実際 − 主部 が深い      : {c4}/{trials}"
          f"  最小の v = {worst['m4']}")
    print(f"    (m5) σ が成分を動かした      : {moved} 件 / 動かさなかった {fixed} 件")
    print(f"    (m5b) v(σf_j − f_j) ≥ 2v(ρ)  : {m5b[0]}/{m5b[1]}"
          f"  最小の v = {m5b[2]}  ← ★Lean の hmv")


def main():
    print("=== ComponentExpansion.lean の主張を塔で測る ===")
    run(3, 3, 60)
    run(3, 4, 20)
    run(5, 3, 12)
    run(7, 2, 12)
    run(2, 4, 60)


if __name__ == "__main__":
    main()
