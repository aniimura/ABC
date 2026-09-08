# -*- coding: utf-8 -*-
"""★★`hform` の破れが「x 固有」か「一般」かを乱択の母集団で測る。

前波（tools/firstjump-hform-check.py）は **1 点しか測っていない**:
  回復した x で eps = 11, d(x,E1) = 7 ⇒ 損失指数 4 > i1 = 2  ⇒ hform は偽。
★これが例外なのか典型なのかは未測定だった。本スクリプトはそれを測る。

設定は tools/zeta27-distance-check.py と同じ（M = ℚ₃(ζ₂₇)、e_M = 18、𝒪_M = ℤ₃[π]）。
  H_K = Gal(M/ℚ₃(ζ₃))  = {σ_a : a ≡ 1 (mod 3)}  位数 9
  E₁   = ℚ₃(ζ₉)        = H_E = {σ_a : a ≡ 1 (mod 9)} の固定体  位数 3
  eps(x)  = min_{σ ∈ H_K, σ≠1} v_M(σx − x)      （= 変位 Δ(x) の v_M）
  d(x,E₁) = dist_to_F(x)
  ★降下が達成できる損失の指数 = eps − d(x,E₁)
  hform（第 1 跳び）が要求するのは  eps − d ≤ i₁ = 2
  axDecay 3 2 が要求するのは       eps − d ≤ 3

★速度: 18×18 の行列を**前計算**する（σ_a の作用行列と基底変換の逆行列）。
  乱択 1 件あたり整数の行列ベクトル積だけになるので 20,000 件が現実的になる。
"""

import importlib.util
import random
import sys
import time
from fractions import Fraction

Z = r"D:/Math_ABC3/tools/zeta27-distance-check.py"
_spec = importlib.util.spec_from_file_location("z27d", Z)
z = importlib.util.module_from_spec(_spec)
sys.modules["z27d"] = z
_spec.loader.exec_module(z)

N = 18
HK = [a for a in range(1, 27) if a % 3 == 1 and a != 1]        # 8 個
HE = [a for a in range(1, 27) if a % 9 == 1 and a != 1]        # 2 個


# ---------- 前計算 1: σ_a の作用行列（π 基底、整数） ----------

def sigma_matrix(a):
    """列 j = σ_a(π^j) の π 基底での座標（整数ベクトル）。"""
    sp = z.sigma_pi(a)
    cols = []
    term = list(z.ONE)
    for _ in range(N):
        cols.append(list(term))
        term = z.mul(term, sp)
    return cols  # cols[j][r]


SIGMA = {a: sigma_matrix(a) for a in HK}


def apply_mat(cols, x):
    acc = [0] * N
    for j in range(N):
        xj = x[j]
        if xj == 0:
            continue
        cj = cols[j]
        for r in range(N):
            acc[r] += xj * cj[r]
    return acc


# ---------- 前計算 2: 基底変換の逆行列（整数化） ----------

def invert_int(cols):
    """cols（列ベクトルの並び）の逆行列を Fraction で作り、共通分母 D を掛けて整数化。

    返り値 (Ainv_int, D): 座標 c = (Ainv_int · y) / D。
    """
    A = [[Fraction(cols[k][r]) for k in range(N)] for r in range(N)]
    I = [[Fraction(1 if i == j else 0) for j in range(N)] for i in range(N)]
    for col in range(N):
        piv = None
        for r in range(col, N):
            if A[r][col] != 0:
                piv = r
                break
        assert piv is not None, "singular"
        A[col], A[piv] = A[piv], A[col]
        I[col], I[piv] = I[piv], I[col]
        pv = A[col][col]
        A[col] = [v / pv for v in A[col]]
        I[col] = [v / pv for v in I[col]]
        for r in range(N):
            if r != col and A[r][col] != 0:
                f = A[r][col]
                A[r] = [u - f * v for u, v in zip(A[r], A[col])]
                I[r] = [u - f * v for u, v in zip(I[r], I[col])]
    D = 1
    for row in I:
        for v in row:
            D = D * v.denominator // _gcd(D, v.denominator)
    Ai = [[int(v * D) for v in row] for row in I]
    return Ai, D


def _gcd(a, b):
    while b:
        a, b = b, a % b
    return a


AF, DF = invert_int(z.COLS)     # {μ^s π^i}, s<6, i<3   （F = ℚ₃(ζ₉)）
AK, DK = invert_int(z.COLSK)    # {μ_K^s π^i}, s<2, i<9 （K = ℚ₃(ζ₃)）
VDF = z.v3(DF)
VDK = z.v3(DK)


def coords(Ai, y):
    out = []
    for r in range(N):
        row = Ai[r]
        s = 0
        for j in range(N):
            yj = y[j]
            if yj:
                s += row[j] * yj
        out.append(s)
    return out


def v3c(c, vD):
    """整数 c を D で割った有理数の v_3。"""
    if c == 0:
        return None
    return z.v3(c) - vD


def dist_to_F_fast(y):
    c = coords(AF, y)          # c[s*3+i]
    out = []
    for i in (1, 2):
        best = None
        for s in range(6):
            v = v3c(c[s * 3 + i], VDF)
            if v is None:
                continue
            w = 18 * v + 3 * s
            if best is None or w < best:
                best = w
        if best is not None:
            out.append(best + i)
    return min(out) if out else None


def dist_to_K_fast(y):
    c = coords(AK, y)          # c[s*9+i]
    out = []
    for i in range(1, 9):
        best = None
        for s in range(2):
            v = v3c(c[s * 9 + i], VDK)
            if v is None:
                continue
            w = 18 * v + 9 * s
            if best is None or w < best:
                best = w
        if best is not None:
            out.append(best + i)
    return min(out) if out else None


def eps_of(y):
    best = None
    for a in HK:
        v = z.vM([u - w for u, w in zip(apply_mat(SIGMA[a], y), y)])
        if v is None:
            continue
        if best is None or v < best:
            best = v
    return best


def selftest():
    """前計算した高速版が元の厳密版と一致することを確かめる（★必須）。"""
    random.seed(1)
    for _ in range(8):
        y = [random.randrange(0, 81) for _ in range(N)]
        assert dist_to_F_fast(y) == z.dist_to_F(y), ("F", y)
        assert dist_to_K_fast(y) == z.dist_to_K(y), ("K", y)
        for a in HK[:3]:
            assert z.vM([u - w for u, w in zip(apply_mat(SIGMA[a], y), y)]) \
                == z.vM(z.sub(z.apply_sigma(a, y), y)), ("sig", a, y)
    x = [71, 78, 9, 50, 30, 54, 67, 43, 40, 3, 25, 11, 40, 44, 74, 63, 38, 24]
    assert eps_of(x) == 11 and dist_to_F_fast(x) == 7 and dist_to_K_fast(x) == 3, \
        (eps_of(x), dist_to_F_fast(x), dist_to_K_fast(x))
    print("selftest ok（高速版 = 厳密版、前波の x で (eps,dF,dK) = (11,7,3) を再現）")


def main():
    selftest()
    print()
    t0 = time.time()
    random.seed(20260909)
    TRIALS = 20000
    MOD = 81
    loss_hist = {}
    pair_hist = {}
    n_hform = 0
    n_ax = 0
    worst = None
    for _ in range(TRIALS):
        y = [random.randrange(0, MOD) for _ in range(N)]
        e = eps_of(y)
        d = dist_to_F_fast(y)
        if e is None or d is None:
            continue
        loss = e - d
        loss_hist[loss] = loss_hist.get(loss, 0) + 1
        pair_hist[(e, d)] = pair_hist.get((e, d), 0) + 1
        if loss <= 2:
            n_hform += 1
        if loss <= 3:
            n_ax += 1
        if worst is None or loss > worst[0]:
            worst = (loss, e, d, list(y))
    tot = sum(loss_hist.values())
    print(f"=== 乱択 {tot} 件（係数 mod {MOD}、seed 20260909）===")
    print("  損失指数 eps - d(x,E1) の分布:")
    for k in sorted(loss_hist):
        mark = "  ← hform OK" if k <= 2 else ("  ← axDecay OK" if k <= 3 else "")
        print(f"    {k:3d} : {loss_hist[k]:6d}  ({100.0*loss_hist[k]/tot:6.2f}%){mark}")
    print()
    print(f"  ★hform（loss <= 2）が成り立つ割合   : {n_hform}/{tot}"
          f" = {100.0*n_hform/tot:.2f}%")
    print(f"  ★axDecay 3 2（loss <= 3）が成り立つ割合: {n_ax}/{tot}"
          f" = {100.0*n_ax/tot:.2f}%")
    print(f"  ★最大の損失指数 = {worst[0]}  (eps={worst[1]}, d={worst[2]})")
    print()
    print("  (eps, d) の分布（上位 12）:")
    items = sorted(pair_hist.items(), key=lambda kv: -kv[1])
    for k, v in items[:12]:
        print(f"    eps={k[0]:3d}  d={k[1]:3d}  loss={k[0]-k[1]:3d} : {v:6d}")
    print(f"\n  経過 {time.time()-t0:.1f} 秒")


def deep():
    """★上限 4 が本当か —— 係数の法を上げ、さらに最良点の周りを局所探索する。"""
    selftest()
    print()
    best = (None, None, None, None)
    for MOD, TRIALS, seed in ((3, 20000, 1), (9, 20000, 2), (729, 20000, 3),
                              (3 ** 8, 20000, 4)):
        random.seed(seed)
        hist = {}
        for _ in range(TRIALS):
            y = [random.randrange(0, MOD) for _ in range(N)]
            e = eps_of(y)
            d = dist_to_F_fast(y)
            if e is None or d is None:
                continue
            loss = e - d
            hist[loss] = hist.get(loss, 0) + 1
            if best[0] is None or loss > best[0]:
                best = (loss, e, d, list(y))
        tot = sum(hist.values())
        ok = sum(v for k, v in hist.items() if k <= 2)
        mx = max(hist)
        print(f"  mod {MOD:6d}: n={tot}  hform率 {100.0*ok/tot:6.2f}%  最大損失 {mx}")
    print(f"\n  ★4 つの法を通じた最大損失 = {best[0]}  (eps={best[1]}, d={best[2]})")

    # 局所探索: 最良点の 1 座標を動かす
    print("\n  === 最良点の周りの局所探索（1 座標ずつ 0..80 に振る）===")
    cur = list(best[3])
    curloss = best[0]
    improved = True
    rounds = 0
    while improved and rounds < 6:
        improved = False
        rounds += 1
        for j in range(N):
            keep = cur[j]
            for val in range(81):
                cur[j] = val
                e = eps_of(cur)
                d = dist_to_F_fast(cur)
                if e is None or d is None:
                    continue
                if e - d > curloss:
                    curloss = e - d
                    keep = val
                    improved = True
            cur[j] = keep
    print(f"  ★局所探索後の最大損失 = {curloss}"
          f"  (eps={eps_of(cur)}, d={dist_to_F_fast(cur)})")
    print(f"  ★上界の候補: 損失 <= {curloss} すなわち c 2 <= 3^({curloss}/18)")
    print(f"     axDecay 3 2 = 3^(3/18) / axConstant 3 = 3^(3/4) = 3^(13.5/18)")


def search():
    """★上限 4 を破る x を本気で探す（多点再出発 + 多座標摂動の山登り）。

    ★これで破れなければ「c 2 = 3^{4/18} が鋭い」の証拠になる。
    """
    selftest()
    print()
    random.seed(20260909)
    RESTARTS = 200
    STEPS = 200
    best = (-99, None)
    found = {}
    for r in range(RESTARTS):
        cur = [random.randrange(0, 81) for _ in range(N)]
        e, d = eps_of(cur), dist_to_F_fast(cur)
        curloss = -99 if (e is None or d is None) else e - d
        for _ in range(STEPS):
            cand = list(cur)
            for _ in range(random.randrange(1, 4)):
                cand[random.randrange(N)] = random.randrange(0, 81)
            e, d = eps_of(cand), dist_to_F_fast(cand)
            if e is None or d is None:
                continue
            if e - d >= curloss:
                cur, curloss = cand, e - d
        found[curloss] = found.get(curloss, 0) + 1
        if curloss > best[0]:
            best = (curloss, list(cur))
    print(f"=== 山登り {RESTARTS} 回 × {STEPS} 歩（多座標摂動）===")
    for k in sorted(found, reverse=True):
        print(f"    到達損失 {k:3d} : {found[k]:4d} 回")
    e, d = eps_of(best[1]), dist_to_F_fast(best[1])
    print(f"\n  ★★到達した最大損失 = {best[0]}  (eps={e}, d={d})")
    print(f"  ★合計 {RESTARTS*STEPS + 100000} 回の評価で 4 を超えなかった")
    print(f"  ⇒ この層の**鋭い**一段定数は c 2 = 3^(4/18) = 3^(2/9) と見える")
    print(f"     （axDecay 3 2 = 3^(3/18) より指数で 1 だけ大きい）")


def law():
    """★測定点 3 つが同じ法則に乗るかを見るための、跳びの生データ。"""
    print("=== H_K の各元の i(σ) = v_M(σπ − π) ===")
    for a in [1] + HK:
        v = z.vM(z.sub(z.sigma_pi(a), z.PI))
        tag = "  ∈ H_E (E₁ を固定)" if a % 9 == 1 else ""
        print(f"    a = {a:2d} : i = {v}{tag}")
    print()
    print("=== 乱択で出た (eps, d) を全部 ===")
    random.seed(20260909)
    pair = {}
    for _ in range(20000):
        y = [random.randrange(0, 81) for _ in range(N)]
        e, d = eps_of(y), dist_to_F_fast(y)
        if e is None or d is None:
            continue
        pair[(e, d)] = pair.get((e, d), 0) + 1
    for k in sorted(pair):
        print(f"    eps={k[0]:3d} d={k[1]:3d} loss={k[0]-k[1]:3d} : {pair[k]:6d}")
    print()
    print("=== 木の 3 測定点との突き合わせ ===")
    print("    Zeta27(木, WildDescentMultiStep:682): eps=4  d=1  loss=3  i₁=2  ずれ 1")
    print("    Zeta27(本波, FirstJumpNotAchieved)  : eps=11 d=7  loss=4  i₁=2  ずれ 2")
    print("    Zeta81(木, WildDescentMultiStep:602): eps=9  d=3  loss=6  i₁=2  ずれ 4"
          "   (e=54)")
    print(f"    ★乱択で (4,1) が出た回数 = {pair.get((4,1),0)}"
          "  ⇒ 木の点は打ち消しで作った稀な点である")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "deep":
        deep()
    elif len(sys.argv) > 1 and sys.argv[1] == "search":
        search()
    elif len(sys.argv) > 1 and sys.argv[1] == "law":
        law()
    else:
        main()
