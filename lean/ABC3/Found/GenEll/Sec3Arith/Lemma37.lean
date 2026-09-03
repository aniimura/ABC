/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Cor43Arith
import ABC3.Meta.Claim
import ABC3.Found.GenEll.Sec3Arith.Lemma35

/-!
# Sec3Arith —— `[GenEll] Lemma 3.7` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll



/-! ## ★★★★★★★★★★★★★★★★`Lemma 3.7`, (a) の数値の核 -/

/-- ★★★★★★★★★★★★★★★★**`Lemma 3.7`, (a) の数値の核**。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

    `100d·(ht^Falt + C·d^ϵ) ≤ l`  ⟹  `d·deg∞ < l·log(2)`

★`Prop 3.4` を `ϵ₀ = 1/6`（`12(1+ϵ₀) = 14`）で使った `deg∞ ≤ 14·ht^Falt + A` と、
`ht^Falt ≥ B`（下に有界）から、**`C ≔ |A| + 100|B| + 1`** と取れば出る。

★★★機構: 求めるのは `A < 55·ht^Falt + 69·C·d^ϵ` である
（`14 d·ht^Falt + dA < 0.69·(100 d·ht^Falt + 100 d C d^ϵ)` を整理した形）。
`ht^Falt ≥ B`・`d^ϵ ≥ 1`・`C = |A| + 100|B| + 1` から
`55·ht^Falt + 69C·d^ϵ ≥ −55|B| + 69|A| + 6900|B| + 69 > |A| ≥ A` ✓。
★`100·log 2 > 69` が効いている。 -/
theorem lemma_3_7_numeric (d degInf htF A B C l eps : ℝ)
    (hd : 1 ≤ d) (heps : 0 < eps)
    (hdeg : degInf ≤ 14 * htF + A) (hB : B ≤ htF)
    (hC : C = |A| + 100 * |B| + 1)
    (hl : 100 * d * (htF + C * d ^ eps) ≤ l)
    (hlog : (0.69 : ℝ) ≤ Real.log 2) :
    d * degInf < l * Real.log 2 := by
  have hdpos : (0:ℝ) < d := by linarith
  have hde : (1:ℝ) ≤ d ^ eps := Real.one_le_rpow hd heps.le
  have hAabs : A ≤ |A| := le_abs_self A
  have hBabs : -|B| ≤ B := neg_abs_le B
  have hCpos : 0 < C := by rw [hC]; positivity
  have hCde : C ≤ C * d ^ eps := by nlinarith
  have hnn : (0:ℝ) ≤ htF + C * d ^ eps := by
    have h1 : B + C ≤ htF + C * d ^ eps := by linarith
    have h2 : (0:ℝ) ≤ B + C := by
      rw [hC]; nlinarith [abs_nonneg A, abs_nonneg B]
    linarith
  have hlpos : (0:ℝ) ≤ l := le_trans (by positivity) hl
  have h1 : d * degInf ≤ d * (14 * htF + A) := mul_le_mul_of_nonneg_left hdeg hdpos.le
  have h2 : (0.69 : ℝ) * l ≤ l * Real.log 2 := by nlinarith
  have h3 : (0.69 : ℝ) * (100 * d * (htF + C * d ^ eps)) ≤ 0.69 * l := by linarith
  have hstrict : A < 55 * htF + 69 * (C * d ^ eps) := by
    have h5 : 69 * C ≤ 69 * (C * d ^ eps) := by nlinarith
    rw [hC] at h5
    nlinarith [abs_nonneg A, abs_nonneg B]
  have hkey : d * A < d * (55 * htF + 69 * (C * d ^ eps)) :=
    mul_lt_mul_of_pos_left hstrict hdpos
  nlinarith [h1, h2, h3, hkey]

def lemma_3_7_numeric.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7, (a)(数値の核——C ≔ |A| + 100|B| + 1)",
    sectionId := "genell-lemma-3-7" }

def lemma_3_7_numeric.needs : List ABC3.Meta.ProofObligation :=
  [ .otherPaper "[GenEll]"
      ("Proposition 3.4(ϵ₀ = 1/6 で deg∞ ≤ 14·ht^Falt + A、および ht^Falt ≥ B)" ++
       "——★第 1・第 2 の合成は §9-980。第 3 の ≲ が残る") 8,
    .otherPaper "[GenEll]"
      "Theorem 3.8(Lemma 3.7 の下流)——★Serre の開像定理。**残る**" 11,
    .implicitStep
      ("★★★★★★測定(2026-08-29): Lemma 3.7, (a) の数値の核は " ++
       "**C ≔ |A| + 100|B| + 1** と取ることに尽きる。" ++
       "★求めるのは A < 55·ht^Falt + 69·C·d^ϵ であり、" ++
       "ht^Falt ≥ B・d^ϵ ≥ 1・C の定義から " ++
       "55·ht^Falt + 69C·d^ϵ ≥ −55|B| + 69|A| + 6900|B| + 69 > |A| ≥ A。" ++
       "★★100·log 2 > 69 が効いている") 6 ]

end ABC3.Found.GenEll
