"""
ABC予想 高速探索エンジン

篩（Sieve）を用いた高速 rad(n) 計算と ABC hit の系統的探索。
quality分布の統計分析、最高quality tripleの発見を行う。
"""

import sys
import io
import math
import time
from collections import defaultdict

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')


def compute_smallest_prime_factor(limit):
    """エラトステネスの篩で各整数の最小素因数を計算する"""
    spf = list(range(limit + 1))  # spf[i] = i の最小素因数
    for i in range(2, int(limit**0.5) + 1):
        if spf[i] == i:  # i は素数
            for j in range(i * i, limit + 1, i):
                if spf[j] == j:
                    spf[j] = i
    return spf


def rad_from_spf(n, spf):
    """最小素因数テーブルを用いた高速 rad(n) 計算"""
    if n <= 1:
        return 1
    result = 1
    temp = n
    while temp > 1:
        p = spf[temp]
        result *= p
        while temp % p == 0:
            temp //= p
    return result


def compute_rad_table(limit, spf):
    """rad(n) のテーブルを一括計算"""
    rad_table = [0] * (limit + 1)
    rad_table[0] = 0
    rad_table[1] = 1
    for n in range(2, limit + 1):
        rad_table[n] = rad_from_spf(n, spf)
    return rad_table


def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def factorize(n, spf):
    """最小素因数テーブルを用いた高速素因数分解"""
    factors = {}
    while n > 1:
        p = spf[n]
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        factors[p] = e
    return factors


def factorization_str(n, spf):
    if n == 1:
        return "1"
    factors = factorize(n, spf)
    parts = []
    for p in sorted(factors):
        if factors[p] == 1:
            parts.append(str(p))
        else:
            parts.append(f"{p}^{factors[p]}")
    return " x ".join(parts)


def search_abc_hits(max_c, spf, rad_table):
    """
    c <= max_c の範囲で ABC hit (quality > 1) を高速探索する。

    最適化: rad(abc) = rad(a) * rad(b) * rad(c) / (共通因子) ではなく、
    gcd(a,b)=1 の条件を使って効率的に計算。
    """
    hits = []
    for c in range(3, max_c + 1):
        rad_c = rad_table[c]
        for a in range(1, c // 2 + 1):
            b = c - a
            if gcd(a, b) != 1:
                continue
            # gcd(a,b) = 1 なので rad(abc) = rad(a)*rad(b)*rad(c) / gcd の重複除去
            # ただし gcd(a,b)=1 => gcd(a,c)=1, gcd(b,c)=1
            # よって rad(abc) の計算は rad(a*b*c) で、共通素因数なし
            rad_abc = rad_table[a] * rad_table[b] * rad_c
            # 但し上の計算は a,b,c が互いに素なので正しい
            if c > rad_abc:  # quality > 1
                q = math.log(c) / math.log(rad_abc)
                hits.append((a, b, c, rad_abc, q))
    hits.sort(key=lambda x: -x[4])
    return hits


def search_high_quality(max_c, spf, rad_table, min_q=0.9):
    """quality >= min_q のABC tripleを探索"""
    results = []
    for c in range(3, max_c + 1):
        rad_c = rad_table[c]
        for a in range(1, c // 2 + 1):
            b = c - a
            if gcd(a, b) != 1:
                continue
            rad_abc = rad_table[a] * rad_table[b] * rad_c
            if rad_abc < 2:
                continue
            q = math.log(c) / math.log(rad_abc)
            if q >= min_q:
                results.append((a, b, c, rad_abc, q))
    results.sort(key=lambda x: -x[4])
    return results


def quality_distribution(max_c, spf, rad_table):
    """quality の分布を計算する"""
    bins = defaultdict(int)
    total = 0
    q_above_1 = 0
    max_q = 0

    for c in range(3, max_c + 1):
        rad_c = rad_table[c]
        for a in range(1, c // 2 + 1):
            b = c - a
            if gcd(a, b) != 1:
                continue
            total += 1
            rad_abc = rad_table[a] * rad_table[b] * rad_c
            if rad_abc < 2:
                continue
            q = math.log(c) / math.log(rad_abc)
            bin_idx = int(q * 10) / 10  # 0.1刻み
            bins[bin_idx] += 1
            if q > 1:
                q_above_1 += 1
            if q > max_q:
                max_q = q

    return bins, total, q_above_1, max_q


def epsilon_analysis(max_c, spf, rad_table):
    """各 epsilon に対して q > 1 + epsilon の triple 数を数える"""
    epsilons = [0.0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5]
    counts = {eps: 0 for eps in epsilons}

    for c in range(3, max_c + 1):
        rad_c = rad_table[c]
        for a in range(1, c // 2 + 1):
            b = c - a
            if gcd(a, b) != 1:
                continue
            rad_abc = rad_table[a] * rad_table[b] * rad_c
            if rad_abc < 2:
                continue
            q = math.log(c) / math.log(rad_abc)
            for eps in epsilons:
                if q > 1 + eps:
                    counts[eps] += 1

    return counts


def main():
    print("=" * 70)
    print("  ABC予想 高速探索エンジン")
    print("  Sieve-based rad computation + systematic ABC hit search")
    print("=" * 70)
    print()

    # Phase 1: 篩の計算
    MAX_C = 10000  # c <= 10^4 (10^5 はメモリ・時間の制約で縮小)
    SIEVE_LIMIT = MAX_C * 3  # a*b*c の最大値に備えて余裕を持つ
    # ただし rad_table は MAX_C まで、a*b*c の rad は個別計算
    # 実際には rad(a)*rad(b)*rad(c) で計算するので rad_table は MAX_C まででOK

    print(f"[Phase 1] 篩の計算 (limit = {MAX_C})...")
    t0 = time.time()
    spf = compute_smallest_prime_factor(MAX_C)
    rad_table = compute_rad_table(MAX_C, spf)
    t1 = time.time()
    print(f"  完了: {t1 - t0:.2f}秒")
    print()

    # Phase 2: ABC hit の探索
    print(f"[Phase 2] ABC hit (quality > 1) の探索 (c <= {MAX_C})...")
    t0 = time.time()
    hits = search_abc_hits(MAX_C, spf, rad_table)
    t1 = time.time()
    print(f"  完了: {t1 - t0:.2f}秒")
    print(f"  ABC hit 総数: {len(hits)}")
    print()

    # 上位20件
    print("  --- 最高 quality ABC triple (上位20件) ---")
    print(f"  {'#':>3} {'a':>10} {'b':>10} {'c':>10} {'rad(abc)':>12} {'quality':>10}")
    print(f"  {'---':>3} {'----------':>10} {'----------':>10} {'----------':>10} {'------------':>12} {'----------':>10}")
    for i, (a, b, c, r, q) in enumerate(hits[:20]):
        print(f"  {i+1:>3} {a:>10} {b:>10} {c:>10} {r:>12} {q:>10.6f}")
    print()

    # 上位5件の詳細分析
    print("  --- 上位5件の素因数分解 ---")
    for i, (a, b, c, r, q) in enumerate(hits[:5]):
        print(f"  #{i+1}: ({a}, {b}, {c}), q = {q:.6f}")
        print(f"    a = {factorization_str(a, spf)}")
        print(f"    b = {factorization_str(b, spf)}")
        print(f"    c = {factorization_str(c, spf)}")
        print(f"    rad(abc) = {r}")
        print()

    # Phase 3: epsilon 解析
    print(f"[Phase 3] epsilon 解析 (c <= {MAX_C})...")
    counts = epsilon_analysis(MAX_C, spf, rad_table)
    print()
    print("  ABC予想の数値的証拠:")
    print("  epsilon が大きいほど、q > 1+eps の triple 数は減少する")
    print()
    print(f"  {'epsilon':>10} {'threshold':>12} {'count':>10} {'ABC予想の予測':>15}")
    print(f"  {'----------':>10} {'------------':>12} {'----------':>10} {'---------------':>15}")
    for eps in sorted(counts.keys()):
        prediction = "有限個" if eps > 0 else "無限個"
        print(f"  {eps:>10.2f} {1+eps:>12.2f} {counts[eps]:>10} {prediction:>15}")

    # Phase 4: quality分布
    print()
    print(f"[Phase 4] quality 分布 (c <= {min(MAX_C, 10000)})...")
    small_limit = min(MAX_C, 10000)
    bins, total, q_above_1, max_q = quality_distribution(small_limit, spf, rad_table)
    print(f"  総 ABC triple 数: {total}")
    print(f"  quality > 1 の数: {q_above_1} ({100*q_above_1/total:.4f}%)")
    print(f"  最大 quality: {max_q:.6f}")
    print()
    print("  quality 分布ヒストグラム:")
    for bin_val in sorted(bins.keys()):
        if bins[bin_val] > 0:
            bar_len = min(int(bins[bin_val] / max(bins.values()) * 50), 50)
            print(f"  [{bin_val:.1f}, {bin_val+0.1:.1f}): {bins[bin_val]:>8} {'#' * bar_len}")

    # Phase 5: ABC hit の c に対する増加率
    print()
    print("[Phase 5] ABC hit 数の増加率分析")
    print()
    thresholds = [100, 500, 1000, 5000, 10000, 50000, MAX_C]
    print(f"  {'c_max':>10} {'hits':>8} {'hits/c_max':>12} {'log比':>10}")
    print(f"  {'----------':>10} {'--------':>8} {'------------':>12} {'----------':>10}")
    for th in thresholds:
        if th > MAX_C:
            break
        n_hits = sum(1 for _, _, c, _, q in hits if c <= th)
        ratio = n_hits / th if th > 0 else 0
        log_ratio = math.log(n_hits) / math.log(th) if n_hits > 0 and th > 1 else 0
        print(f"  {th:>10} {n_hits:>8} {ratio:>12.6f} {log_ratio:>10.4f}")

    print()
    print("=" * 70)
    print("  分析完了")
    print("=" * 70)


if __name__ == "__main__":
    main()
