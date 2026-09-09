# -*- coding: utf-8 -*-
"""★`GainedTowerDescent.lean:94` が「証明していない」と書いた閉じた形を測る。

同ファイルの docstring（測定 4）:
  ★閉じた形も出た（★4,194 本 → 不一致 0 件だが本ファイルは証明していない）:
     Λ_k = t_k + Σ_{j<k} t_j − t_1·(p^{k−1}−1)/(p−1)

漸化式（同 :184）:
  Λ_0 = 0
  Λ_{m+1} = max( A_{m+1} , max(0, A_{m+1} − t_1) + p·Λ_m ),  A_{m+1} = t_{m+1} − (p−1)·J_m
  J_m = t_1 + ⋯ + t_m

★本スクリプトは (i) 閉じた形が**一般には偽**であることを示し、
(ii) **どの条件の下で真か**を突き止める。
"""

import itertools
import sys


def jumpSum(t, m):
    return sum(t[j] for j in range(1, m + 1))


def gainedLoss(p, t, k):
    L = 0
    for m in range(1, k + 1):
        A = t[m] - (p - 1) * jumpSum(t, m - 1)
        L = max(A, max(0, A - t[1]) + p * L)
    return L


def closed(p, t, k):
    g = sum(p ** i for i in range(k - 1))          # (p^{k-1}-1)/(p-1)
    return t[k] + jumpSum(t, k - 1) - t[1] * g


def branch_flags(p, t, k):
    """各段で (第2枝が勝ったか, max(0,·) が正だったか) を返す。"""
    L = 0
    out = []
    for m in range(1, k + 1):
        A = t[m] - (p - 1) * jumpSum(t, m - 1)
        b2 = max(0, A - t[1]) + p * L
        out.append((b2 >= A, A - t[1] >= 0))
        L = max(A, b2)
    return out


def scan(p, kmax, vals):
    tot = same = 0
    diff_examples = []
    cond_ok_and_same = cond_ok_and_diff = 0
    for k in range(1, kmax + 1):
        for tail in itertools.product(vals, repeat=k):
            t = [0] + list(tail)
            if any(t[j] <= 0 for j in range(1, k + 1)):
                continue
            tot += 1
            L = gainedLoss(p, t, k)
            C = closed(p, t, k)
            fl = branch_flags(p, t, k)
            # ★候補の条件: 2 段目以降で「第 2 枝が勝ち、かつ A_m ≥ t_1」
            cond = all(b2 and pos for (b2, pos) in fl[1:])
            if L == C:
                same += 1
                if cond:
                    cond_ok_and_same += 1
            else:
                if cond:
                    cond_ok_and_diff += 1
                if len(diff_examples) < 4:
                    diff_examples.append((k, tail, L, C, fl))
    return tot, same, diff_examples, cond_ok_and_same, cond_ok_and_diff


def main():
    vals = list(range(1, 13))
    print("=== 閉じた形 Λ_k = t_k + J_{k−1} − t_1·(p^{k−1}−1)/(p−1) の検査 ===")
    for p in (2, 3, 5):
        tot, same, ex, cok, cdiff = scan(p, 4, vals)
        print(f"  p={p}: 全 {tot} 本中 一致 {same} 本"
              f"（{100.0*same/tot:.1f}%）  ★不一致 {tot-same} 本")
        print(f"     ★候補条件（2 段目以降で第 2 枝かつ A_m ≥ t_1）を満たすもの: "
              f"一致 {cok} / 不一致 {cdiff}")
        for (k, tail, L, C, fl) in ex[:2]:
            print(f"     不一致例 k={k} t={tail}  Λ={L}  閉じた形={C}  枝={fl}")


if __name__ == "__main__":
    main()
