/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.PrimesOfSize
import ABC3.Found.GenEll.Elementary
import ABC3.Meta.Claim
import ABC3.Found.GenEll.Cor43Arith.Corollary43

/-!
# Cor43Arith —— `[GenEll] Corollary 4.4` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll



/-! ## ★★★★★★★★★★★★`Corollary 4.4`（`h = 0`） -/

/-- ★★★★★★★★★★★★**`Corollary 4.4, (c)` の `l∘` の評価**。

原文 (GenEll p.23):
> Corollary 4.4.

    `l∘ ≤ 23040 · 100d · ht^Falt([E_L]) + 2x_S + C · d`

★原文は『`The proof is entirely similar to [but slightly easier than] the proof of
Corollary 4.3, except that ... when applying Lemma 4.1, we take “h” to be 0`』と書く。
★★`h = 0` なので `l ≤ 2·x_A` だけになり、係数は `2·3·12 = 72 ≤ 100` で済む
——`d^{1+ϵ}` ではなく **`d`** で足りるのはこのためである。 -/
theorem cor_4_4_arith (d xS xSo dgInf htF B C'' : ℝ)
    (hd : 1 ≤ d) (hC'' : 0 ≤ C'') (hB : B ≤ htF)
    (hxSo : xSo ≤ xS + (5 / 2) * (23040 * d) * dgInf)
    (hdeg : dgInf ≤ (12 * (6 / 5)) * htF + C'')
    (l : ℝ) (hl : l ≤ 2 * xSo) :
    l ≤ 23040 * 100 * d * htF + 2 * xS + (23040 * (5 * C'' + 28 * |B|)) * d := by
  have hdpos : (0:ℝ) < d := by linarith
  have hdegd : (23040 * d) * dgInf ≤ (23040 * d) * ((12 * (6 / 5)) * htF + C'') :=
    mul_le_mul_of_nonneg_left hdeg (by positivity)
  have hBd : 28 * (23040 * d) * B ≤ 28 * (23040 * d) * htF :=
    mul_le_mul_of_nonneg_left hB (by positivity)
  have hBabs : -(28 * 23040 * |B|) * d ≤ 28 * (23040 * d) * B := by
    have h1 : -|B| ≤ B := neg_abs_le B
    nlinarith [abs_nonneg B]
  nlinarith [hl, hxSo, hdegd, hBd, hBabs]

/-- ★★★★★★★★★★**`Corollary 4.4, (c)` の `l•` の評価**（`log-diff` つき）。

原文 (GenEll p.23):
> Corollary 4.4.

    `l• ≤ 23040 · 100d · ht^Falt + 6d · log-diff_Mell + 2x_S + C · d` -/
theorem cor_4_4_arith_bullet (d xS xSb ldiff htF B C'' : ℝ)
    (hd : 1 ≤ d) (hC'' : 0 ≤ C'') (hB : B ≤ htF)
    (hxSb : xSb ≤ xS + (3 * 12 * 23040) * d * htF + 3 * d * ldiff + d * C'')
    (l : ℝ) (hl : l ≤ 2 * xSb) :
    l ≤ 23040 * 100 * d * htF + 6 * d * ldiff + 2 * xS
        + (23040 * (2 * C'' + 28 * |B|)) * d := by
  have hdpos : (0:ℝ) < d := by linarith
  have hBd : 28 * (23040 * d) * B ≤ 28 * (23040 * d) * htF :=
    mul_le_mul_of_nonneg_left hB (by positivity)
  have hBabs : -(28 * 23040 * |B|) * d ≤ 28 * (23040 * d) * B := by
    have h1 : -|B| ≤ B := neg_abs_le B
    nlinarith [abs_nonneg B]
  nlinarith [hl, hxSb, hBd, hBabs]

def cor_4_4_arith.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 23,
    item := "Corollary 4.4, (c)(l∘ の評価——h = 0 で d で足りる)",
    sectionId := "genell-cor-4-4" }

def cor_4_4_arith_bullet.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 23,
    item := "Corollary 4.4, (c)(l• の評価——log-diff つき)",
    sectionId := "genell-cor-4-4" }

def cor_4_4_arith.needs : List ABC3.Meta.ProofObligation :=
  [ .otherPaper "[GenEll]"
      "Theorem 3.8, (b)(Galois 表現が全射)——★Serre の開像定理。**残る**" 11,
    .implicitStep
      ("★★★★★測定(2026-08-29): 原文は Corollary 4.4 を『entirely similar to " ++
       "[but slightly easier than] the proof of Corollary 4.3』とし、" ++
       "違いは Theorem 3.8 の条件 (a) ではなく (b) を使うことと、" ++
       "Lemma 4.1 で h = 0 と取ることだけである。" ++
       "★★h = 0 なので l ≤ 2·x_A だけになり、係数は 2·3·12 = 72 ≤ 100 で済む" ++
       "——**d^{1+ϵ} ではなく d で足りる**のはこのためである") 6 ]

end ABC3.Found.GenEll
