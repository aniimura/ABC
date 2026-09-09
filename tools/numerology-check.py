# -*- coding: utf-8 -*-
"""★`LossExponentMatch.lean` の数合わせを塔で測る（規律 9: 式には母集団）。

測るのは 7 点:

  (n1) v(w) = p か（Lean は `‖B‖ = ‖π‖^{p + v(f_{j₀})}` と書いている）
  (n2) p ∣ v(f_j) か（`𝒪_{E₁}` の付値は p の倍数。Lean の `hdvd`）
  (n3) d = min_{1≤j≤p−1}( v(f_j) + j ) が `Layer.dist_to_E1` と一致するか
  (n4) ★j₀ ≠ jstar は起きるか（★起きるなら `ResidueSeparation.loss_le_of_residue_ne` の
       等号仮定 `‖B‖ = ‖π‖^{d+(p−1)}` は具体層で**満たされない**。一般化が要る根拠）
  (n5) v(B) = v(w) + v(f_{j₀}) か（先頭の対が打ち消さない、の付値版）
  (n6) ★v(R) % p ≠ j₀ % p か（Lean の `hv`/`hne`。★残る唯一の入力）
  (n7) ★v(B) + j₀ ≤ d + 2p−2 と、実際の v(σx−x) ≤ d + 2p−2（＝ loss ≤ 2p−2）
"""

import importlib.util
import random
import sys

K3 = r"D:/Math_ABC3/tools/hform-k3-check.py"
_spec = importlib.util.spec_from_file_location("hfk3num", K3)
hf = importlib.util.module_from_spec(_spec)
sys.modules["hfk3num"] = hf
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


def smul(k, a):
    return [k * x for x in a]


def run(p, n, trials, seed=20260909):
    L = hf.Layer(p, n)
    F = L.F
    vD = L.N * L.vD                          # v_L(D)
    iv = {a: F.v(F.sub(F.sigma_pi(a), F.PI)) for a in L.HK}
    a0 = [a for a in L.HK if iv[a] == L.i1 + 1][0]
    b = (a0 - 1) // p
    one_plus_pi = F.red([1, 1] + [0] * (F.N - 2))
    w = F.sub(F.powm(one_plus_pi, p * b), F.ONE)
    vw = F.v(w)

    random.seed(seed)
    n2 = [0, 0]
    n3 = [0, 0]
    n4 = 0                                    # j₀ ≠ jstar の件数
    n5 = [0, 0]
    n6 = [0, 0]
    n7a = [0, 0]
    n7b = [0, 0]
    n8 = [0, 0, 0, 0]
    n10 = [0, 0]
    n11 = [0, 0, None]           # [v(f_j) >= 0 の件数, 総数, 最小の v(f_j)]
    n9 = [0, 0, None]            # [成立, 総数, 最大の loss]            # [誤差が深い, ★同位, 誤差が浅い, 総数]
    worst = None
    for _ in range(trials):
        y = [p ** random.randrange(0, 4) * random.randrange(0, p ** 3)
             for _ in range(L.N)]
        cx = L.coords(y)
        Df = [elt_from_coords(L, cx, i) for i in range(p)]
        Df.append([0] * L.N)
        vf = []
        for i in range(p + 1):
            z = F.v(Df[i])
            vf.append(None if z is None else z - vD)

        tail = [j for j in range(1, p) if vf[j] is not None]
        if not tail:
            continue
        m = min(vf[j] for j in tail)
        j0 = max(j for j in tail if vf[j] == m)
        jstar = min(tail, key=lambda j: vf[j] + j)
        d = vf[jstar] + jstar

        # (n11) ★係数の整数性: x が整なら v(f_j) ≥ 0 か
        #      （CoefficientIntegrality.norm_coeff_le_one_of_valK が言っていること）
        vx = F.v(y)
        if vx is not None and vx >= 0:
            for j in range(p):
                if vf[j] is None:
                    continue
                n11[1] += 1
                n11[0] += (vf[j] >= 0)
                if n11[2] is None or vf[j] < n11[2]:
                    n11[2] = vf[j]
        # (n2)
        for j in tail:
            n2[1] += 1
            n2[0] += (vf[j] % p == 0)
        # (n3)
        d_ref = L.dist_to_E1(y)
        n3[1] += 1
        n3[0] += (d_ref == d)
        # (n4)
        if j0 != jstar:
            n4 += 1
        # (n5)
        Belt = F.mul(w, add(smul(j0, Df[j0]), smul(j0 + 1, Df[j0 + 1])))
        vB = F.v(Belt)
        vB = None if vB is None else vB - vD
        n5[1] += 1
        n5[0] += (vB is not None and vB == vw + vf[j0])
        # (n6) R = w Σ_{j ≠ j0} ( j f_j + (j+1) f_{j+1} ) π^j
        Relt = [0] * L.N
        for j in range(p):
            if j == j0:
                continue
            coef = add(smul(j, Df[j]), smul(j + 1, Df[j + 1]))
            Relt = add(Relt, F.mul(coef, F.powm(F.PI, j)))
        Relt = F.mul(w, Relt)
        vR = F.v(Relt)
        vR = None if vR is None else vR - vD
        n6[1] += 1
        n6[0] += (vR is None or vR % p != j0 % p)
        # (n7)
        if vB is not None:
            n7a[1] += 1
            n7a[0] += (vB + j0 <= d + 2 * p - 2)
        # (n8) 誤差項 err = (σx − x) − (B π^{j₀} + R) の位を主部と比べる
        main = add(F.mul(Belt, F.powm(F.PI, j0)), Relt)
        actual = [0] * L.N
        for j in range(p):
            actual = add(actual, F.mul(L._apply(L.sig[a0], Df[j]),
                                       F.powm(F.sigma_pi(a0), j)))
            actual = [u - z for u, z in zip(actual, F.mul(Df[j], F.powm(F.PI, j)))]
        err = [u - z for u, z in zip(actual, main)]
        vmain = F.v(main)
        verr = F.v(err)
        n8[3] += 1
        if verr is None or (vmain is not None and verr > vmain):
            n8[0] += 1                        # 誤差の方が深い（strict 版が使える）
        elif vmain is not None and verr == vmain:
            n8[1] += 1                        # ★同じ位（ここだけが危ない）
        else:
            n8[2] += 1                        # 誤差の方が浅い（結論は出る）

        # (n10) 誤差が「係数の大きさ × ‖ρ‖²」以下か: v(err) ≥ v(f_{j₀}) + 2p
        n10[1] += 1
        if verr is None or verr - vD >= vf[j0] + 2 * p:
            n10[0] += 1

        veps = F.v([u - z for u, z in zip(L._apply(L.sig[a0], y), y)])
        if veps is not None:
            n7b[1] += 1
            n7b[0] += (veps <= d + 2 * p - 2)
            loss = veps - d
            if worst is None or loss > worst:
                worst = loss
        # (n9) ★真の loss（σ を H_K 全体で動かした min）
        tl = L.loss(y)
        if tl is not None:
            n9[1] += 1
            n9[0] += (tl <= 2 * p - 2)
            if n9[2] is None or tl > n9[2]:
                n9[2] = tl

    print(f"  p={p} n={n} 標本 {trials}  v(w)={vw}  p={p}  2p−2={2 * p - 2}")
    print(f"    (n1) v(w) = p                : {'成立' if vw == p else '★不成立'}")
    print(f"    (n11) ★x 整 ⇒ v(f_j) ≥ 0     : {n11[0]}/{n11[1]}"
          f"  最小の v(f_j) = {n11[2]}")
    print(f"    (n2) p | v(f_j)              : {n2[0]}/{n2[1]}")
    print(f"    (n3) d の 2 通りの計算が一致 : {n3[0]}/{n3[1]}")
    print(f"    (n4) ★j₀ ≠ jstar             : {n4}/{n3[1]} 件"
          f"  ← 等号仮定が破れる標本")
    print(f"    (n5) v(B) = v(w) + v(f_{{j₀}}) : {n5[0]}/{n5[1]}")
    print(f"    (n6) ★v(R) % p ≠ j₀ % p      : {n6[0]}/{n6[1]}  ← Lean の hv/hne")
    print(f"    (n7) v(B) + j₀ ≤ d + 2p−2    : {n7a[0]}/{n7a[1]}")
    print(f"         実際の v(σx−x) ≤ d+2p−2 : {n7b[0]}/{n7b[1]}"
          f"  ★測った loss の最大 = {worst}")
    print(f"    (n8) 誤差 err と主部の位     : 深い {n8[0]} / ★同位 {n8[1]}"
          f" / 浅い {n8[2]}  （総数 {n8[3]}）")
    print(f"    (n10) v(err) ≥ v(f_j₀)+2p    : {n10[0]}/{n10[1]}"
          f"  ← loss_le_of_error_small の仮定")
    print(f"    (n9) ★真の loss ≤ 2p−2       : {n9[0]}/{n9[1]}"
          f"  最大の loss = {n9[2]}  ← σ を H_K 全体で動かした min")


def main():
    print("=== LossExponentMatch.lean の数合わせを塔で測る ===")
    run(3, 3, 200)
    run(3, 4, 60)
    run(5, 3, 40)
    run(7, 2, 30)
    run(2, 4, 200)


if __name__ == "__main__":
    main()
