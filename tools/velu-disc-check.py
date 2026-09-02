# -*- coding: utf-8 -*-
"""Velu の商の判別式の恒等式を数値で確かめる(2026-09-02、第 1384)。

    Delta(E)^l = Delta(E/C) * ( prod_{P in C\{O}} (2 y_P + a1 x_P + a3) )^4

l = 3, 5, 7 について Tate 標準形の族で厳密に確かめる。
使い方:  'C:/Users/Aruta/miniforge3/envs/py311env/python.exe' tools/velu-disc-check.py
"""
from fractions import Fraction as F


def b_invs(a1, a2, a3, a4, a6):
    b2 = a1 * a1 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3 * a3 + 4 * a6
    b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4
    return b2, b4, b6, b8


def disc(a1, a2, a3, a4, a6):
    b2, b4, b6, b8 = b_invs(a1, a2, a3, a4, a6)
    return -b2 * b2 * b8 - 8 * b4 ** 3 - 27 * b6 * b6 + 9 * b2 * b4 * b6


def add(P, Q, a1, a2, a3, a4, a6):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2 and y1 + y2 + a1 * x2 + a3 == 0:
        return None
    if P == Q:
        lam = (3 * x1 * x1 + 2 * a2 * x1 + a4 - a1 * y1) / (2 * y1 + a1 * x1 + a3)
    else:
        lam = (y2 - y1) / (x2 - x1)
    nu = y1 - lam * x1
    x3 = lam * lam + a1 * lam - a2 - x1 - x2
    y3 = -(lam + a1) * x3 - nu - a3
    return (x3, y3)


def velu(a1, a2, a3, a4, a6, kernel):
    reps = []
    seen = set()
    for (x, y) in kernel:
        negy = -y - a1 * x - a3
        if (x, negy) in seen:
            continue
        seen.add((x, y))
        reps.append((x, y))
    v = F(0)
    w = F(0)
    b2, _, _, _ = b_invs(a1, a2, a3, a4, a6)
    prod = F(1)
    for (x, y) in kernel:
        prod *= (2 * y + a1 * x + a3)
    for (x, y) in reps:
        gx = 3 * x * x + 2 * a2 * x + a4 - a1 * y
        gy = -2 * y - a1 * x - a3
        vP = gx if (2 * y + a1 * x + a3) == 0 else 2 * gx - a1 * gy
        v += vP
        w += gy * gy + x * vP
    return (a1, a2, a3, a4 - 5 * v, a6 - b2 * v - 7 * w), prod


def test(name, a1, a2, a3, a4, a6, P, l):
    D = disc(a1, a2, a3, a4, a6)
    pts = []
    Q = P
    for _ in range(1, l):
        pts.append(Q)
        Q = add(Q, P, a1, a2, a3, a4, a6)
    if Q is not None:
        print(name, "order != ", l)
        return False
    (c1, c2, c3, c4p, c6p), prod = velu(a1, a2, a3, a4, a6, pts)
    Dp = disc(c1, c2, c3, c4p, c6p)
    ok = (D ** l == Dp * prod ** 4)
    print("%-22s l=%d  identity=%s" % (name, l, ok))
    return ok


def main():
    allok = True
    # l = 3 : y^2 + a1 x y + a3 y = x^3, P = (0,0)
    for a1v, a3v in [(F(1), F(1)), (F(2), F(1)), (F(1), F(3)), (F(-2), F(5))]:
        allok &= test("l3 a1=%s a3=%s" % (a1v, a3v), a1v, F(0), a3v, F(0), F(0),
                      (F(0), F(0)), 3)
    # l = 5 : y^2 + (1-c) x y - c y = x^3 - c x^2, P = (0,0)
    for c in [F(1), F(2), F(3), F(-1), F(5)]:
        allok &= test("l5 c=%s" % c, 1 - c, -c, -c, F(0), F(0), (F(0), F(0)), 5)
    # l = 7 : b = d^3 - d^2, c = d^2 - d
    for d in [F(2), F(3), F(-1), F(5)]:
        b = d ** 3 - d ** 2
        c = d ** 2 - d
        allok &= test("l7 d=%s" % d, 1 - c, -b, -b, F(0), F(0), (F(0), F(0)), 7)
    print("ALL OK" if allok else "FAILED")


if __name__ == "__main__":
    main()
