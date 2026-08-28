/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Cor43Arith
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★`Lemma 3.5`／`Lemma 3.7` の数値の核（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17–p.19。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★これは何か

§3 の `Lemma 3.5`／`Lemma 3.7` は、`Interface/GenEll/EllModuli.lean` の
`EllModuliData`（witness の無い界面）の上に立っている。
★★**その中の「数値の核」だけは界面に依らずに取れる**——本ファイルはそれを取る。

### `Lemma 3.5`（原文 p.17）

    `(1/(12(1+ϵ)))·l·deg∞(E) ≤ ht^Falt(E) + 2log(l) + C`

★原文の証明は 3 つの入力の**割り算**である:
`deg∞(E_H) = l·deg∞(E)`（`Lemma 3.2, (ii)`）、
`deg∞ ≤ 12(1+ϵ)·ht^Falt + K`（`Prop 3.4`）、
`ht^Falt(E_H) ≤ ht^Falt(E) + 2log(l) + C₀`（[FC] Chapter I, Prop 2.7 と `(1,1)`-形式の積分）。

### `Lemma 3.7`, (a)（原文 p.18）

    `100d·(ht^Falt + C·d^ϵ) ≤ l`  ⟹  `d·deg∞ < l·log(2)`

★`Prop 3.4` を `ϵ₀ = 1/6`（すなわち `12(1+ϵ₀) = 14`）で使って `deg∞ ≤ 14·ht^Falt + A`、
`ht^Falt ≥ B`（下に有界）とし、**`C ≔ |A| + 100|B| + 1`** と取ればよい。
★★`100·log 2 > 69` と `d^ϵ ≥ 1` が効く。

## ★★★何が残るか（測定、2026-08-29）

| 入力 | 状態 |
|---|---|
| ★数値の核（本ファイル） | ★★**取れた** |
| `Prop 3.4`（`deg∞ ≲ 12(1+ϵ)ht^Falt`） | △ 第 1・第 2 の合成は `§9-980`、第 3 が残る |
| `Lemma 3.2, (ii)`（`deg∞(E_H) = l·deg∞(E)`） | ☆残る（Tate 曲線） |
| `ht^Falt(E_H) ≤ ht^Falt(E) + 2log(l)` | ☆残る（[FC] Chapter I, Prop 2.7） |
| `Theorem 3.8`（`Lemma 3.7` の下流） | ☆残る（Serre の開像定理） |

★`.src` は条つき——指標には数えない。
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

/-! ## ★出典の紐付け(`.src`)——★**条つきである。指標には数えない** -/

def lemma_3_5_numeric.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(数値の核——3 つの入力の割り算)",
    sectionId := "genell-lemma-3-5" }

def lemma_3_7_numeric.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7, (a)(数値の核——C ≔ |A| + 100|B| + 1)",
    sectionId := "genell-lemma-3-7" }

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
