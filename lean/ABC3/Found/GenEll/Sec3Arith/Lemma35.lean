/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Cor43Arith
import ABC3.Meta.Claim

/-!
# Sec3Arith —— `[GenEll] Lemma 3.5` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll



/-! ## ★★★★★★★★★★★★`Lemma 3.5` の数値の核 -/

/-- ★★★★★★★★★★★★**`Lemma 3.5` の数値の核**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

    `(1/(12(1+ϵ)))·l·deg∞ ≤ ht^Falt + 2log(l) + C₀ + K/(12(1+ϵ))`

★原文の証明は 3 つの入力の**割り算**である——
`deg∞(E_H) = l·deg∞(E)`・`deg∞ ≤ 12(1+ϵ)ht^Falt + K`・
`ht^Falt(E_H) ≤ ht^Falt(E) + 2log(l) + C₀`。

★★**定数が `ϵ` に依り `E, F, H_F, l` に依らない**ことは、
`K`・`C₀` が `l` を含まないことに現れている
——★原文が明記している量化子の順序である。 -/
theorem lemma_3_5_numeric (eps l degInf degInf' htF htF' K C0 : ℝ)
    (heps : 0 < eps)
    (hdeg' : degInf' = l * degInf)
    (hchain : degInf' ≤ 12 * (1 + eps) * htF' + K)
    (hfal : htF' ≤ htF + 2 * Real.log l + C0) :
    (1 / (12 * (1 + eps))) * l * degInf
      ≤ htF + 2 * Real.log l + C0 + K / (12 * (1 + eps)) := by
  have hpos : (0:ℝ) < 12 * (1 + eps) := by nlinarith
  have h1 : l * degInf ≤ 12 * (1 + eps) * (htF + 2 * Real.log l + C0) + K := by
    rw [← hdeg']
    nlinarith [hchain, hfal, hpos]
  rw [mul_assoc, one_div, inv_mul_eq_div, div_le_iff₀ hpos]
  field_simp
  nlinarith [h1]

/-! ## ★出典の紐付け(`.src`)——★**条つきである。指標には数えない** -/

def lemma_3_5_numeric.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(数値の核——3 つの入力の割り算)",
    sectionId := "genell-lemma-3-5" }

def lemma_3_5_numeric.needs : List ABC3.Meta.ProofObligation :=
  [ .otherPaper "[GenEll]"
      ("Lemma 3.2, (ii)(deg∞(E_H) = l·deg∞(E))——★Tate 曲線。**残る**") 8,
    .otherPaper "[GenEll]"
      ("Proposition 3.4(deg∞ ≲ 12(1+ϵ)ht^Falt)" ++
       "——★第 1・第 2 の合成は §9-980。第 3 の ≲ が残る") 8,
    .otherPaper "[FC]"
      ("Chapter I, Proposition 2.7(等因子写像が延びること)と (1,1)-形式の積分" ++
       "——ht^Falt(E_H) ≤ ht^Falt(E) + 2log(l)。**残る**") 6,
    .implicitStep
      ("★★★★★測定(2026-08-29): Lemma 3.5 の証明は 3 つの入力の**割り算**である。" ++
       "★本ファイルはその数値の核を界面(EllModuliData)に依らずに取った。" ++
       "★★定数が ϵ に依り E, F, H_F, l に依らないことは、K・C₀ が l を含まないことに現れている" ++
       "——原文が明記している量化子の順序である") 5 ]

end ABC3.Found.GenEll
