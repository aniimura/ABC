# -*- coding: utf-8 -*-
"""K0 の最適値を厳密整数で測る。

模型: 円分塔 L = Q_p(zeta_{p^k}), E1 = Q_p(zeta_{p^{k-1}})
  t k = p            (第 1 跳び)
  m k = p^(k-1)
  e k = p^(k-2) * (p-1)   (= e_{E1}, phi(p^{k-1}))
必要な不等式 (LossToExit / FirstJumpRoute の hjump):
  (p-1) * (t + (p-2)) <= e
Nat の切り捨て引き算に合わせて max(x,0) で書く。
"""

def nsub(a, b):
    return a - b if a > b else 0

def npow(p, n):
    return p ** n

def needed(p, k):
    """(p-1)*(t+(p-2)) <= p^(k-2)*(p-1) が成り立つか"""
    lhs = nsub(p, 1) * (p + nsub(p, 2))
    rhs = npow(p, nsub(k, 2)) * nsub(p, 1)
    return lhs <= rhs, lhs, rhs

primes = [2, 3, 5, 7, 11, 13, 101]
print("p   最小の k で成立   (k=0..3 の可否)")
for p in primes:
    firsts = None
    row = []
    for k in range(0, 9):
        ok, lhs, rhs = needed(p, k)
        if k <= 4:
            row.append("k=%d:%s(%d<=%d)" % (k, "T" if ok else "F", lhs, rhs))
        if ok and firsts is None:
            # 以後ずっと成立するか確認
            if all(needed(p, kk)[0] for kk in range(k, 40)):
                firsts = k
    print("p=%-4d 最小 K0 = %s" % (p, firsts))
    print("      " + "  ".join(row))

print()
print("--- 2*(p-1) <= p^(k-2) の形（抽象核）---")
for p in primes:
    firsts = None
    for k in range(0, 40):
        if 2 * nsub(p, 1) <= npow(p, nsub(k, 2)):
            if all(2 * nsub(p, 1) <= npow(p, nsub(kk, 2)) for kk in range(k, 40)):
                firsts = k
                break
    print("p=%-4d 最小 K0 = %s" % (p, firsts))

print()
print("--- p=2 で余白が消えるか: (p-1)*(t+(p-2)) vs (p-1)*t ---")
for p in [2, 3, 5]:
    for t in [1, 2, 3, 5]:
        a = nsub(p, 1) * (t + nsub(p, 2))
        b = nsub(p, 1) * t
        print("p=%d t=%d : slackあり %d  / 木の形 %d  / 差 %d" % (p, t, a, b, a - b))
