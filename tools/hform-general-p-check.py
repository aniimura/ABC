# -*- coding: utf-8 -*-
"""★測った最大損失 4 が「2·i₁」なのか「2(p−1)」なのか「p に依らない 4」なのかを決める。

前波（tools/hform-k3-check.py）は p = 3 の 2 つの深さしか測っていない:
  k = 2 (e=18) も k = 3 (e=54) も  i₁ = 2、axDecay = 3 目盛り、★最大損失 = 4 目盛り。
★p = 3 では p − 1 = 2 = i₁ なので、4 = 2·i₁ = 2(p−1) = 4 の 3 通りが区別できない。

⇒ ★p = 2 と p = 5 で測れば 1 回で決まる。

  p = 2: p − 1 = 1
  p = 5: p − 1 = 4
"""

import importlib.util
import random
import sys
import time

K3 = r"D:/Math_ABC3/tools/hform-k3-check.py"
_spec = importlib.util.spec_from_file_location("hfk3", K3)
hf = importlib.util.module_from_spec(_spec)
sys.modules["hfk3"] = hf
_spec.loader.exec_module(hf)


def measure(p, n, trials, seed=20260909, mod=None):
    t0 = time.time()
    L = hf.Layer(p, n)
    mod = mod or p ** 4
    random.seed(seed)
    hist = {}
    worst = None
    for _ in range(trials):
        y = [random.randrange(0, mod) for _ in range(L.N)]
        l = L.loss(y)
        if l is None:
            continue
        hist[l] = hist.get(l, 0) + 1
        if worst is None or l > worst[0]:
            worst = (l, L.eps(y), L.dist_to_E1(y))
    tot = sum(hist.values())
    ok = sum(v for k, v in hist.items() if k <= L.i1)
    mx = max(hist)
    gen = max(hist, key=lambda k: hist[k])
    # axDecay p (n-1) は層の付値でちょうど p 目盛り（HformPopulation の抽象核）
    print(f"--- p={p} n={n}  e={L.N}  |H_K|={len(L.HK)+1}  i1={L.i1}  p-1={p-1} ---")
    print(f"    axDecay(目盛り) = {p}   生成的な損失 = {gen}   ★最大損失 = {mx}"
          f"   (eps={worst[1]}, d={worst[2]})")
    print(f"    hform率 = {100.0*ok/tot:.2f}%   n={tot}   分布 = "
          f"{dict(sorted(hist.items()))}")
    print(f"    ★比: 最大/i1 = {mx}/{L.i1}   最大/(p-1) = {mx}/{p-1}"
          f"   最大-axDecay = {mx-p}   経過 {time.time()-t0:.1f}s")
    return dict(p=p, n=n, e=L.N, i1=L.i1, gen=gen, mx=mx)


def main():
    rows = []
    # p = 3 の再現（前波と一致するはず）
    rows.append(measure(3, 2, 20000))
    rows.append(measure(3, 3, 20000))
    # p = 2
    rows.append(measure(2, 3, 20000))
    rows.append(measure(2, 4, 20000))
    rows.append(measure(2, 5, 8000))
    # p = 5
    rows.append(measure(5, 2, 8000))
    print()
    print("=== まとめ ===")
    print("  p  n   e   i1  p-1  生成的  最大  最大=2·i1?  最大=2(p-1)?  最大=4?")
    for r in rows:
        p, i1 = r["p"], r["i1"]
        print(f"  {p}  {r['n']}  {r['e']:3d}  {i1:2d}  {p-1:3d}"
              f"  {r['gen']:5d}  {r['mx']:4d}"
              f"   {str(r['mx'] == 2*i1):10s} {str(r['mx'] == 2*(p-1)):12s}"
              f" {str(r['mx'] == 4)}")


def climb(p, n, restarts, steps, seed=20260909, mod=None):
    """★稀な最大値を山登りで探す（乱択だけでは p が大きいと届かない）。"""
    t0 = time.time()
    L = hf.Layer(p, n)
    print(f"  Layer({p},{n}) 構築 {time.time()-t0:.1f}s   e={L.N}  i1={L.i1}"
          f"  |H_K|={len(L.HK)+1}")
    mod = mod or p ** 4
    random.seed(seed)
    best = None
    reach = {}
    for _ in range(restarts):
        cur = [random.randrange(0, mod) for _ in range(L.N)]
        cl = L.loss(cur)
        if cl is None:
            cl = -10 ** 6
        for _ in range(steps):
            cand = list(cur)
            for _ in range(random.randrange(1, 4)):
                cand[random.randrange(L.N)] = random.randrange(0, mod)
            l = L.loss(cand)
            if l is not None and l >= cl:
                cur, cl = cand, l
        reach[cl] = reach.get(cl, 0) + 1
        if best is None or cl > best:
            best = cl
    print(f"  到達分布 {dict(sorted(reach.items(), reverse=True))}")
    print(f"  ★★到達した最大損失 = {best}   予測 p+1 = {p+1}"
          f"   一致 {best == p+1}   経過 {time.time()-t0:.1f}s")
    return best


def biased(p, n, trials, seed=20260909):
    """★係数の p 進付値を散らして測る（一様乱択だと p が大きいと全部 loss = i₁ になる）。

    y[j] = p^{r_j} · u  （r_j は 0..4 の一様、u は 0..p^4 の一様）。
    """
    t0 = time.time()
    L = hf.Layer(p, n)
    random.seed(seed)
    hist = {}
    worst = None
    for _ in range(trials):
        y = [p ** random.randrange(0, 5) * random.randrange(0, p ** 4)
             for _ in range(L.N)]
        l = L.loss(y)
        if l is None:
            continue
        hist[l] = hist.get(l, 0) + 1
        if worst is None or l > worst[0]:
            worst = (l, L.eps(y), L.dist_to_E1(y))
    tot = sum(hist.values())
    ok = sum(v for k, v in hist.items() if k <= L.i1)
    print(f"--- 偏り付き  p={p} n={n}  e={L.N}  i1={L.i1}  axDecay={p} ---")
    print(f"    分布 = {dict(sorted(hist.items()))}")
    print(f"    hform率 = {100.0*ok/tot:.2f}%   ★最大損失 = {max(hist)}"
          f"   (eps={worst[1]}, d={worst[2]})   予測 p+1 = {p+1}")
    print(f"    経過 {time.time()-t0:.1f}s")
    return max(hist)


def predict():
    """★法則 max = p + 1 の予言を p = 5 で試す（p=2,3 では既に成立）。"""
    print("=== 予言の試験: max loss = p + 1 ===")
    print("[対照] p=3 n=3 （既知の答え 4）")
    climb(3, 3, 60, 120)
    print("[対照] p=2 n=4 （既知の答え 3）")
    climb(2, 4, 60, 120)
    print("[★本番] p=5 n=3 （予言 6）")
    climb(5, 3, int(sys.argv[2]) if len(sys.argv) > 2 else 40,
          int(sys.argv[3]) if len(sys.argv) > 3 else 80)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "predict":
        predict()
    else:
        main()
