# -*- coding: utf-8 -*-
"""Gained ルートの 1 層の値段は「単項 axDecay p k」に落ちるかを厳密整数で測る。

GainedTowerDescent.lean:
  gainedLoss p t 0       = 0
  gainedLoss p t (m+1)   = max(A_{m+1}, max(0, A_{m+1} - t 1) + p * gainedLoss p t m)
    ただし A_{m+1} = t(m+1) - (p-1) * jumpSum t m,  jumpSum t m = t1 + ... + tm

rpow_le_prod_axDecay (:458) + towerBudget_iff (:435):
  p^{J/(p^k e)} <= ∏_{i in Icc 1 k} axDecay p i   <=>  (p-1)^2 * J <= (p^{k+1} - p) * e

単項版(自分で導出、下でも数値検算する):
  axDecay p k = p^{(1/(p-1)) p^{1-k}}
  p^{J/(p^k e)} <= axDecay p k  <=>  J/(p^k e) <= p^{1-k}/(p-1)
                                <=>  (p-1) * J <= p * e
"""
from fractions import Fraction

def jump_sum(t, m):
    return sum(t[j] for j in range(1, m + 1))

def gained_loss(p, t, m):
    lam = 0
    for i in range(1, m + 1):
        A = t[i] - (p - 1) * jump_sum(t, i - 1)
        lam = max(A, max(0, A - t[1]) + p * lam)
    return lam

def prod_ok(p, e, k, J):
    return (p - 1) ** 2 * J <= (p ** (k + 1) - p) * e

def single_ok(p, e, J):
    return (p - 1) * J <= p * e

# 木が §10 ChainCheck で挙げている実データ
# (name, p, e(=底の絶対分岐指数), t[1..n])
CASES = [
    ("Q3(zeta27)/Q3(zeta3)  n=2", 3, 2, {1: 2, 2: 8}),
    ("Q3(zeta81)/Q3(zeta9)  n=2", 3, 6, {1: 8, 2: 26}),
    ("Q3(zeta81)/Q3(zeta3)  n=3", 3, 2, {1: 2, 2: 8, 3: 26}),
    ("p=5 cyclotomic       n=2", 5, 4, {1: 4, 2: 24}),
    ("p=5 cyclotomic       n=3", 5, 4, {1: 4, 2: 24, 3: 124}),
    ("p=2 cyclotomic       n=2", 2, 1, {1: 1, 2: 3}),
    ("p=2 cyclotomic       n=3", 2, 1, {1: 1, 2: 3, 3: 7}),
]

print("case                          n   Lambda   積の条件           単項の条件")
print("                                          (p-1)^2 L <= (p^{n+1}-p)e   (p-1)L <= p e")
for name, p, e, t in CASES:
    n = max(t)
    L = gained_loss(p, t, n)
    lhs_p, rhs_p = (p - 1) ** 2 * L, (p ** (n + 1) - p) * e
    lhs_s, rhs_s = (p - 1) * L, p * e
    print("%-28s %d  %6d   %5d <= %-5d  %s   %5d <= %-4d %s" %
          (name, n, L, lhs_p, rhs_p, "OK " if lhs_p <= rhs_p else "NG ",
           lhs_s, rhs_s, "OK" if lhs_s <= rhs_s else "★NG"))

print()
print("--- 2 つの閾値の比 ---")
print("積の条件が許す J の上限 : (p^{n+1}-p) e / (p-1)^2")
print("単項が許す J の上限     : p e / (p-1)")
print("比 = (p^{n+1}-p)/((p-1) p) = (p^n - 1)/(p-1) = 1 + p + ... + p^{n-1}")
for p in (2, 3, 5):
    for n in (1, 2, 3, 4):
        r = Fraction(p ** n - 1, p - 1)
        print("  p=%d n=%d : 比 = %s" % (p, n, r))

print()
print("--- 単項で押さえるのに必要な指数 (実データ) ---")
print("必要:  J/(p^n e) <= X  ,  axDecay p n の指数は (1/(p-1)) p^{1-n}")
for name, p, e, t in CASES:
    n = max(t)
    L = gained_loss(p, t, n)
    need = Fraction(L, p ** n * e)
    have = Fraction(1, p - 1) * Fraction(1, p ** (n - 1))
    print("%-28s 必要 %-10s  axDecay p %d の指数 %-10s  %s" %
          (name, str(need), n, str(have), "OK" if need <= have else "★足りない"))
