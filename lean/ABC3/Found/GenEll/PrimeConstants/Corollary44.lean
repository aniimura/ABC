/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.Tactic.NormNum.Prime
import ABC3.Found.GenEll.PrimeNumberTheorem
import ABC3.Found.GenEll.PrimeConstants.Corollary43

/-!
# PrimeConstants —— `[GenEll] Corollary 4.4` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll

open Real Finset

/-- ★★★★**`Corollary 4.4` の (c)**——`Lemma 4.1` の `h` を `0` に取る場合。

原文 (GenEll p.23):
> Corollary 4.3, except that instead of applying condition (a) of Theorem 3.8, we

★原文:『Also, when applying Lemma 4.1, we take “`h`” to be **`0`**』。
★★`8h` の項が消えるので係数は `2·3·12 = 72 ≤ **100**` で足り、
末項も `C·d^{1+ϵ}` ではなく **`C·d`** になる——これが `Corollary 4.3` との差である。 -/
theorem cor44_numeric (d F dinf xS xbad xT extra l A B P : ℝ)
    (hd1 : 1 ≤ d) (hP1 : 1 ≤ P) (hdP : d ≤ P)
    (hA : dinf ≤ 12 * (1 + 1/5) * F + A)
    (hxbad : xbad ≤ 5/2 * (23040 * (d * dinf)))
    (hB : -B ≤ F) (hBnn : 0 ≤ B) (hxT : 0 ≤ xT)
    (hl : l ≤ 2 * (xS + xbad + xT + extra)) :
    l ≤ 23040 * 100 * (d * F) + 2 * extra + 2 * xS
        + (23040 * (5 * |A| + 28 * B) + 2 * xT + 1) * P := by
  have hdpos : (0:ℝ) < d := by linarith
  have hAabs : A ≤ |A| := le_abs_self A
  have hdA : d * dinf ≤ 12 * (1 + 1/5) * (d * F) + d * A := by
    have hm := mul_le_mul_of_nonneg_left hA hdpos.le
    linarith [hm]
  have h1 : d * A ≤ d * |A| := mul_le_mul_of_nonneg_left hAabs hdpos.le
  have h2 : d * |A| ≤ P * |A| := mul_le_mul_of_nonneg_right hdP (abs_nonneg A)
  have hdBP : d * B ≤ P * B := mul_le_mul_of_nonneg_right hdP hBnn
  have hxTP : xT ≤ xT * P := by nlinarith
  have hFB : (0:ℝ) ≤ d * F + d * B := by nlinarith
  have hPnn : (0:ℝ) ≤ P := by linarith
  linarith [hl, hxbad, hdA, h1, h2, hdBP, hxTP, hFB, hPnn]

def cor44_numeric.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 23,
    item := "Corollary 4.4 の証明(when applying Lemma 4.1, we take “h” to be 0)",
    sectionId := "genell-cor-4-4" }

end ABC3.Found.GenEll
