"""
フレイ曲線解析スクリプト

ABC triple (a, b, c) に対応するフレイ曲線 E: y^2 = x(x-a)(x+b) の
判別式、導手、Szpiro比を計算し、ABC予想との関係を数値的に調査する。
"""

import sys
import io
import math
from collections import defaultdict

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')


def prime_factors(n):
    """n の素因数とその指数を返す"""
    factors = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def rad(n):
    """n の根基 rad(n) = 素因数の積"""
    if n <= 1:
        return 1
    result = 1
    for p in prime_factors(n):
        result *= p
    return result


def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def ord_p(n, p):
    """p 進付値 v_p(n): n が p で割れる回数"""
    if n == 0:
        return float('inf')
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def factorization_str(n):
    """素因数分解の文字列表現"""
    if n == 1:
        return "1"
    factors = prime_factors(n)
    parts = []
    for p in sorted(factors):
        if factors[p] == 1:
            parts.append(str(p))
        else:
            parts.append(f"{p}^{factors[p]}")
    return " x ".join(parts)


# ============================================================
# フレイ曲線の構成と解析
# ============================================================

def frey_curve_discriminant(a, b, c):
    """
    フレイ曲線 E_{a,b,c}: y^2 = x(x - a)(x + b) の判別式を計算

    a + b = c, gcd(a, b) = 1 を仮定。
    判別式: Delta = 16 * (abc)^2 / (gcd factors)
    最小判別式: Delta_min = (abc)^2 / 256 (正規化後)

    実際には:
    E: y^2 = x(x - a)(x + b) の判別式は
    Delta = 16 * a^2 * b^2 * c^2
    (ただし a + b = c の場合)
    """
    # y^2 = x(x-a)(x+b) = x^3 + (b-a)x^2 - ab*x
    # 一般ワイエルシュトラス形式 y^2 = x^3 + Ax^2 + Bx + C
    # A = b - a, B = -a*b, C = 0
    # 判別式 = -4A^3C + A^2B^2 + 18ABC - 4B^3 - 27C^2
    # C = 0 なので: Delta_formula = A^2 * B^2 - 4 * B^3
    # = B^2 * (A^2 - 4B) = (ab)^2 * ((b-a)^2 + 4ab) = (ab)^2 * (a+b)^2 = (abc)^2

    # ただし標準的な定義では係数16がつく:
    # Delta = 16 * (a*b*c)^2  (非最小モデルの場合)

    return 16 * (a * b * c) ** 2


def frey_curve_conductor_factors(a, b, c):
    """
    フレイ曲線の導手 N の素因数分解（簡略版）

    導手 N = prod_{p | abc} f_p  ここで f_p は局所導手指数。

    gcd(a,b) = 1 かつ a+b = c のとき：
    - p が奇素数: f_p = p (半安定還元、導手指数 = 1)
    - p = 2: f_p の計算はより複雑（2での還元型に依存）

    簡略化: N = rad(abc) (半安定の場合の近似)
    より正確には: N は rad(abc) を割り切り、rad(abc) は N を割り切る（2の寄与を除く）
    """
    # 簡略版: 全ての素因数で半安定と仮定
    # 実際のフレイ曲線では Ribet の定理等により
    # 導手は rad(abc) の約数（2を除く）
    abc = a * b * c
    factors = prime_factors(abc)

    conductor_factors = {}
    for p, e in factors.items():
        if p == 2:
            # 2 での導手指数は複雑
            # 大まかに: a,b,c の 2-adic 構造に依存
            # 簡略化として指数を推定
            conductor_factors[p] = min(e, 6)  # 上界
        else:
            # 奇素数: 半安定還元なので導手指数 = 1
            conductor_factors[p] = 1

    return conductor_factors


def compute_conductor(conductor_factors):
    """導手の値を計算"""
    N = 1
    for p, e in conductor_factors.items():
        N *= p ** e
    return N


def szpiro_ratio(a, b, c):
    """
    Szpiro比: log|Delta_min| / log(N)

    Szpiro予想: 任意の eps > 0 に対し、
    log|Delta_min| <= (6 + eps) * log(N) + C(eps)

    ここでは簡略的に:
    - Delta ~ (abc)^2
    - N ~ rad(abc)
    として Szpiro比 = log((abc)^2) / log(rad(abc)) を計算
    """
    abc = a * b * c
    r = rad(abc)
    if r <= 1:
        return 0
    return math.log(abc ** 2) / math.log(r)


def quality(a, b, c):
    """ABC quality: log(c) / log(rad(abc))"""
    r = rad(a * b * c)
    if r <= 1:
        return 0
    return math.log(c) / math.log(r)


# ============================================================
# メイン解析
# ============================================================

def find_abc_triples(max_c, min_quality=1.0):
    """quality >= min_quality の ABC triple を探索"""
    triples = []
    for c in range(3, max_c + 1):
        for a in range(1, c // 2 + 1):
            b = c - a
            if gcd(a, b) != 1:
                continue
            r = rad(a * b * c)
            if r < 2:
                continue
            q = math.log(c) / math.log(r)
            if q >= min_quality:
                triples.append((a, b, c, q))
    triples.sort(key=lambda x: -x[3])
    return triples


def main():
    print("=" * 70)
    print("  フレイ曲線解析 - ABC triple と楕円曲線の関係")
    print("=" * 70)
    print()

    # Phase 1: 有名な ABC triple のフレイ曲線解析
    print("[Phase 1] 有名な ABC triple のフレイ曲線解析")
    print()

    famous_triples = [
        (1, 8, 9),        # 1 + 2^3 = 3^2
        (5, 27, 32),      # 5 + 3^3 = 2^5
        (1, 80, 81),      # 1 + 2^4 * 5 = 3^4
        (32, 49, 81),     # 2^5 + 7^2 = 3^4
        (1, 2, 3),        # 最小の ABC triple (q < 1)
        (1, 48, 49),      # 1 + 2^4 * 3 = 7^2
        (1, 63, 64),      # 1 + 3^2 * 7 = 2^6
        (1, 4374, 4375),  # 高い quality
        (2, 6436341, 6436343),  # 非常に高い quality
    ]

    print(f"  {'(a, b, c)':>30} | {'q':>8} | {'Delta':>18} | {'N(approx)':>12} | {'Szpiro':>8}")
    print(f"  {'-' * 30}-+-{'-' * 8}-+-{'-' * 18}-+-{'-' * 12}-+-{'-' * 8}")

    for a, b, c in famous_triples:
        if a + b != c:
            continue
        if gcd(a, b) != 1:
            continue
        q = quality(a, b, c)
        delta = frey_curve_discriminant(a, b, c)
        cond_factors = frey_curve_conductor_factors(a, b, c)
        N = compute_conductor(cond_factors)
        sz = szpiro_ratio(a, b, c)
        print(f"  ({a}, {b}, {c}){' ' * max(0, 25 - len(f'({a}, {b}, {c})'))} | {q:>8.4f} | {delta:>18} | {N:>12} | {sz:>8.4f}")

    print()

    # Phase 2: 各 ABC triple の詳細解析
    print("[Phase 2] 高 quality ABC triple の詳細フレイ曲線解析")
    print()

    detail_triples = [
        (1, 8, 9),
        (5, 27, 32),
        (1, 80, 81),
        (1, 4374, 4375),
    ]

    for a, b, c in detail_triples:
        print(f"  --- ABC triple: ({a}, {b}, {c}) ---")
        print(f"    a = {factorization_str(a)}")
        print(f"    b = {factorization_str(b)}")
        print(f"    c = {factorization_str(c)}")
        print(f"    a + b = c: {a + b == c}")
        print(f"    gcd(a, b) = {gcd(a, b)}")
        print()

        q = quality(a, b, c)
        r = rad(a * b * c)
        print(f"    rad(abc) = rad({a * b * c}) = {r}")
        print(f"    ABC quality = log({c}) / log({r}) = {q:.6f}")
        print()

        # フレイ曲線
        print(f"    フレイ曲線: y^2 = x(x - {a})(x + {b})")
        print(f"    = y^2 = x^3 + {b - a}*x^2 + {-a * b}*x")
        delta = frey_curve_discriminant(a, b, c)
        print(f"    判別式 Delta = 16 * ({a} * {b} * {c})^2 = {delta}")
        print(f"    log|Delta| = {math.log(abs(delta)):.4f}")
        print()

        # 導手
        cond_factors = frey_curve_conductor_factors(a, b, c)
        N = compute_conductor(cond_factors)
        print(f"    導手の素因数: {cond_factors}")
        print(f"    導手 N (近似) = {N}")
        print(f"    rad(abc) = {r}")
        print(f"    N / rad(abc) = {N / r:.4f}" if r > 0 else "    N / rad(abc) = N/A")
        print(f"    log(N) = {math.log(N):.4f}")
        print()

        # Szpiro 比
        sz = szpiro_ratio(a, b, c)
        print(f"    Szpiro比 = log(Delta) / log(N) = {math.log(delta):.4f} / {math.log(N):.4f} = {sz:.4f}")
        print(f"    Szpiro予想の上界: 6 + eps")
        print(f"    {'[OK] Szpiro比 < 6' if sz < 6 else '[!] Szpiro比 >= 6 (注目)'}")
        print()

        # quality と Szpiro比の関係
        # Szpiro比 ~ 2 * quality (大まかな関係)
        print(f"    quality q = {q:.6f}")
        print(f"    Szpiro比 / (2 * q) = {sz / (2 * q):.6f}" if q > 0 else "")
        print(f"    (理論的には Szpiro比 ~ 2 * quality + 定数)")
        print()

    # Phase 3: Szpiro 比と quality の相関分析
    print("[Phase 3] Szpiro比と quality の相関分析 (c <= 10000)")
    print()

    triples = find_abc_triples(10000, min_quality=0.8)
    print(f"  quality >= 0.8 の triple 数: {len(triples)}")

    # quality vs Szpiro ratio の散布データ
    q_values = []
    sz_values = []
    for a, b, c, q in triples:
        sz = szpiro_ratio(a, b, c)
        q_values.append(q)
        sz_values.append(sz)

    if q_values:
        # 線形回帰 (最小二乗法)
        n = len(q_values)
        sum_q = sum(q_values)
        sum_sz = sum(sz_values)
        sum_q2 = sum(q ** 2 for q in q_values)
        sum_qsz = sum(q * sz for q, sz in zip(q_values, sz_values))

        denom = n * sum_q2 - sum_q ** 2
        if abs(denom) > 1e-10:
            slope = (n * sum_qsz - sum_q * sum_sz) / denom
            intercept = (sum_sz - slope * sum_q) / n
            print(f"  線形回帰: Szpiro比 = {slope:.4f} * quality + {intercept:.4f}")
            print(f"  (理論予測: Szpiro比 ~ 2 * quality)")
            print(f"  実測された傾き / 2 = {slope / 2:.4f}")
        print()

        # quality > 1 の triple について Szpiro比の分布
        abc_hits = [(a, b, c, q) for a, b, c, q in triples if q > 1.0]
        print(f"  ABC hit (quality > 1) の数: {len(abc_hits)}")
        if abc_hits:
            print()
            print(f"  {'(a, b, c)':>25} | {'quality':>8} | {'Szpiro':>8} | Szpiro/6")
            print(f"  {'-' * 25}-+-{'-' * 8}-+-{'-' * 8}-+--------")
            for a, b, c, q in abc_hits[:15]:
                sz = szpiro_ratio(a, b, c)
                print(f"  ({a}, {b}, {c}){' ' * max(0, 20 - len(f'({a}, {b}, {c})'))} | {q:>8.4f} | {sz:>8.4f} | {sz / 6:>6.4f}")

    # Phase 4: j-不変量の解析
    print()
    print("[Phase 4] フレイ曲線の j-不変量")
    print()
    print("  j-不変量 j(E) = 1728 * (4A^3) / Delta  (短ワイエルシュトラス形式)")
    print("  フレイ曲線 y^2 = x(x-a)(x+b) では:")
    print("  完全な j-不変量の計算には変数変換が必要。")
    print("  ここでは c4, c6 を用いた標準的な計算を行う。")
    print()

    for a, b, c in [(1, 8, 9), (5, 27, 32), (1, 80, 81)]:
        # y^2 = x^3 + (b-a)x^2 - ab*x
        # a2 = b - a, a4 = -ab, a1 = a3 = a6 = 0
        a2_coeff = b - a
        a4_coeff = -a * b

        # b2 = 4*a2, b4 = 2*a4
        b2 = 4 * a2_coeff
        b4 = 2 * a4_coeff

        # c4 = b2^2 - 24*b4
        c4 = b2 ** 2 - 24 * b4

        # Delta_disc (from Silverman's formula)
        # c6 = -b2^3 + 36*b2*b4
        c6 = -b2 ** 3 + 36 * b2 * b4

        delta_check = (c4 ** 3 - c6 ** 2) // 1728

        j_inv = c4 ** 3 / delta_check if delta_check != 0 else float('inf')

        print(f"  ({a}, {b}, {c}):")
        print(f"    a2 = {a2_coeff}, a4 = {a4_coeff}")
        print(f"    c4 = {c4}, c6 = {c6}")
        print(f"    Delta = {delta_check}")
        print(f"    j = c4^3 / Delta = {c4 ** 3} / {delta_check} = {j_inv:.4f}")
        print(f"    j-不変量は楕円曲線の同型類を決定する")
        print()

    # Phase 5: ABC予想と Szpiro予想の等価性の数値的検証
    print("[Phase 5] ABC予想 <=> Szpiro予想 の等価性の数値的検証")
    print()
    print("  定理: ABC予想は Szpiro予想と等価である。")
    print("  具体的には:")
    print("  ABC予想: c < K_eps * rad(abc)^{1+eps}")
    print("  Szpiro予想: |Delta_min| < C_eps * N^{6+eps}")
    print()
    print("  質的関係: quality q と Szpiro比 sigma の間に")
    print("  sigma ~ 6q (大まかな関係)")
    print()

    if q_values and sz_values:
        # 相関係数の計算
        mean_q = sum(q_values) / len(q_values)
        mean_sz = sum(sz_values) / len(sz_values)
        cov = sum((q - mean_q) * (sz - mean_sz) for q, sz in zip(q_values, sz_values)) / len(q_values)
        var_q = sum((q - mean_q) ** 2 for q in q_values) / len(q_values)
        var_sz = sum((sz - mean_sz) ** 2 for sz in sz_values) / len(sz_values)

        if var_q > 0 and var_sz > 0:
            corr = cov / (var_q ** 0.5 * var_sz ** 0.5)
            print(f"  quality と Szpiro比の相関係数: {corr:.6f}")
            print(f"  (1に近いほど強い正の相関)")
        print()

        # 各 quality 範囲での Szpiro比の平均
        ranges = [(0.8, 0.9), (0.9, 1.0), (1.0, 1.1), (1.1, 1.2), (1.2, 1.5)]
        print(f"  {'quality 範囲':>15} | {'count':>6} | {'Szpiro比の平均':>14} | {'Szpiro / (6*q_mid)':>18}")
        print(f"  {'-' * 15}-+-{'-' * 6}-+-{'-' * 14}-+-{'-' * 18}")
        for lo, hi in ranges:
            subset = [(q, sz) for q, sz in zip(q_values, sz_values) if lo <= q < hi]
            if subset:
                avg_sz = sum(sz for _, sz in subset) / len(subset)
                q_mid = (lo + hi) / 2
                print(f"  [{lo:.1f}, {hi:.1f}){' ' * 5} | {len(subset):>6} | {avg_sz:>14.4f} | {avg_sz / (6 * q_mid):>18.4f}")

    print()
    print("=" * 70)
    print("  フレイ曲線解析完了")
    print("=" * 70)


if __name__ == "__main__":
    main()
