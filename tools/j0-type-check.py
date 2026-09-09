# -*- coding: utf-8 -*-
"""★Lean に書いた `j₀` の型が、機械の言うことと一致するかを測る。

Lean 側（`lean/ABC3/Found/PGC/MaxMinIndex.lean`）の型はこう言っている:

  `j₀` := `[1, p-1]` のうち `v(f_j)` を最小にする**最大の**添字。このとき
    (1) `∀ j ∈ [1,p-1], ‖f j‖ ≤ ‖f j₀‖`
    (2) `‖ j₀·f_{j₀} + (j₀+1)·f_{j₀+1} ‖ = ‖f_{j₀}‖`（先頭の対が打ち消さない）
    仮定は `p ∤ j₀`（これは `1 ≤ j₀ ≤ p-1` から出る）と `f_p = 0` と
    「尾に非零成分が 1 つある」だけ。

本スクリプトが測るのは 4 点:
  (A) 尾 `[1,p-1]` が丸ごと零になる標本はどれだけあるか（＝型の仮定 `hane` が
      いつ効かないか。★これを測らずに「仮定は軽い」と書かない）
  (B) 尾に零成分が混じる標本はどれだけあるか（★付値版の `hv` が排除している集合）
  (C) `v(f_{j₀+1}) > v(f_{j₀})` が実際に成り立つか（型の定義通りのはず。定義の
      実装可能性の確認）
  (D) ★`v(P_{j₀}) = v(w) + v(f_{j₀})`（＝ Lean の結論 (2) を付値に訳したもの）と
      `v(B_{j₀} − P_{j₀}) > v(P_{j₀})`（＝成分の公式が `j₀` で真）
"""

import importlib.util
import random
import sys

K3 = r"D:/Math_ABC3/tools/hform-k3-check.py"
_spec = importlib.util.spec_from_file_location("hfk3j0", K3)
hf = importlib.util.module_from_spec(_spec)
sys.modules["hfk3j0"] = hf
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


def add(a, b):
    return [x + y for x, y in zip(a, b)]


def smul(k, a):
    return [k * x for x in a]


def run(p, n, trials, seed=20260909):
    L = hf.Layer(p, n)
    F = L.F
    iv = {a: F.v(F.sub(F.sigma_pi(a), F.PI)) for a in L.HK}
    firstj = [a for a in L.HK if iv[a] == L.i1 + 1]
    random.seed(seed)
    a0 = firstj[0]
    b = (a0 - 1) // p
    one_plus_pi = F.red([1, 1] + [0] * (F.N - 2))
    w = F.sub(F.powm(one_plus_pi, p * b), F.ONE)
    vw = F.v(w)

    tail_all_zero = 0        # (A)
    tail_some_zero = 0       # (B)
    c_ok = c_bad = 0         # (C)
    d_pair_ok = d_pair_bad = 0   # (D) 先頭の対が打ち消さないか
    d_form_ok = d_form_bad = 0   # (D) 成分の公式が j₀ で真か
    ex = []
    for _ in range(trials):
        y = [p ** random.randrange(0, 4) * random.randrange(0, p ** 3)
             for _ in range(L.N)]
        cx = L.coords(y)
        Df = [elt_from_coords(L, cx, i) for i in range(p)]
        Df.append([0] * L.N)                       # f_p = 0
        vf = [F.v(Df[i]) for i in range(p + 1)]    # None = 零 = ∞

        tail = [j for j in range(1, p) if vf[j] is not None]
        if not tail:
            tail_all_zero += 1
            continue
        if len(tail) < p - 1:
            tail_some_zero += 1

        m = min(vf[j] for j in tail)
        j0 = max(j for j in tail if vf[j] == m)     # ★Lean の定義そのもの

        # (C) v(f_{j0+1}) > v(f_{j0})  （f_p = 0 は ∞ とみなす）
        vnext = vf[j0 + 1]
        if vnext is None or vnext > vf[j0]:
            c_ok += 1
        else:
            c_bad += 1
            if len(ex) < 3:
                ex.append(("C", j0, vf[j0], vnext))

        # (D) 先頭の対 P = w( j0 f_{j0} + (j0+1) f_{j0+1} )
        inner = add(smul(j0, Df[j0]), smul(j0 + 1, Df[j0 + 1]))
        vinner = F.v(inner)
        if vinner is not None and vinner == vf[j0]:
            d_pair_ok += 1
        else:
            d_pair_bad += 1
            if len(ex) < 6:
                ex.append(("Dpair", j0, vf[j0], vinner))

        DP = F.mul(w, inner)
        vP = F.v(DP)
        diff = [u - v for u, v in zip(L._apply(L.sig[a0], y), y)]
        cb = L.coords(diff)
        DB = elt_from_coords(L, cb, j0)
        if vP is None:
            d_form_bad += 1
        else:
            vR = F.v([u - v for u, v in zip(DB, DP)])
            if vR is None or vR > vP:
                d_form_ok += 1
            else:
                d_form_bad += 1
                if len(ex) < 9:
                    ex.append(("Dform", j0, vP, vR))

    tot = trials
    print(f"  p={p} n={n} 標本 {tot} 件  (v(w)={vw})")
    print(f"    (A) 尾が丸ごと零        : {tail_all_zero} 件")
    print(f"    (B) 尾に零成分が混じる  : {tail_some_zero} 件"
          f"  ← 付値版 hv が排除している集合")
    print(f"    (C) v(f_j0+1) > v(f_j0) : 成立 {c_ok} / 不成立 {c_bad}")
    print(f"    (D) 先頭の対が非打消し  : 成立 {d_pair_ok} / 不成立 {d_pair_bad}")
    print(f"    (D) 成分の公式が j0 で真: 成立 {d_form_ok} / 不成立 {d_form_bad}")
    if ex:
        print(f"    反例: {ex}")


def designed(p, n, trials=30, seed=20260909):
    """★成分が零になる `y` を**設計して**作り、2 つの型の差が本物かを測る。

    座標 `c` の「成分 jz のブロック」を丸ごと 0 にして `y = Σ c_k · basis_k` を作る。
    このとき `f_jz = 0` なので ★付値版（`hv : ‖f j‖ = ‖π‖^{v j}`）は**適用できない**。
    ノルム版（`exists_no_cancel_index`）は適用できる。差が空でないことを示す。
    """
    L = hf.Layer(p, n)
    F = L.F
    cols = L._basis_cols()
    random.seed(seed)
    made = ok = allzero = 0
    for _ in range(trials):
        jz = random.randrange(1, p)          # 零にする尾の添字
        c = [random.randrange(0, p ** 3) for _ in range(L.N)]
        for s in range(L.S):
            c[s * p + jz] = 0
        y = [0] * L.N
        for k in range(L.N):
            if c[k]:
                y = add(y, smul(c[k], cols[k]))
        cx = L.coords(y)
        Df = [elt_from_coords(L, cx, i) for i in range(p)]
        Df.append([0] * L.N)
        vf = [F.v(Df[i]) for i in range(p + 1)]
        if vf[jz] is not None:
            continue                          # 設計に失敗（起きないはず）
        made += 1
        tail = [j for j in range(1, p) if vf[j] is not None]
        if not tail:
            allzero += 1
            continue
        m = min(vf[j] for j in tail)
        j0 = max(j for j in tail if vf[j] == m)
        inner = add(smul(j0, Df[j0]), smul(j0 + 1, Df[j0 + 1]))
        vinner = F.v(inner)
        if vinner is not None and vinner == vf[j0]:
            ok += 1
    print(f"  p={p} n={n} ★設計: f_jz = 0 を作れた {made}/{trials} 件"
          f"（尾が丸ごと零 {allzero} 件）"
          f"  そのうち先頭の対が非打消し {ok} 件")


def main():
    print("=== Lean の j0 の型 vs 機械 ===")
    run(3, 3, 200)
    run(3, 4, 60)
    run(5, 3, 40)
    run(7, 2, 40)
    run(2, 4, 200)
    print("=== 零成分を設計して作る（付値版とノルム版の差） ===")
    designed(3, 3)
    designed(5, 3, 20)


if __name__ == "__main__":
    main()
