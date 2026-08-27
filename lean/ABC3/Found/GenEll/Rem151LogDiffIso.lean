/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.MinFieldCovering
import ABC3.Found.GenEll.LogDiffCongr

/-!
# [GenEll] Remark 1.5.1 の前半 —— **`log-diff` は同型で不変**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.8。

原文 (GenEll p.8):
> function log-diffX on X(Q) [cf. Definition 1.5, (iii)] depends only on the scheme

## ★★★`Remark 1.5.1` は 2 つの主張である——混同しない

| 主張 | 何が要るか | 状態 |
|---|---|---|
| (A) `log-diff_X` は `X_ℚ` にしか依らない | ★**同型不変性だけ**——spreading out は要らない | ★**本ファイル** |
| (B) `log-cond_D` の BD-類は `(X_ℚ, D_ℚ)` にしか依らない | 算術的 spreading out ＋ Σ 上の一様上界 | 計算部分は `LogCondSigma.lean` の `abs_logCond_sub_le_sum_log`。★spreading out は未 |

★★**(A) に spreading out が要らないのは、`log-diff` が点の最小定義体だけで決まるから**である
（`Definition 1.5, (iii)`）。★因子も計量も見ないので、ℤ-モデルの取り方が入り込む余地がない。

## ★★★★証明は `minField_comp_le` を両向きに当てるだけ

`MinFieldCovering.lean` の

    minField F (x_F ≫ φ) ≤ minField F x_F   （任意の射 φ）

を `φ = e.hom` と `φ = e.inv` に当てると、同型の場合は**等号**になる。
★★あとは `logDiffOfField_congr`（`LogDiffCongr.lean`：同型な数体は同じ `log-diff`）を継ぐ。

## ★★★★★対比——一般の射では等号にならない

`minField_comp_le` は**不等号**であり、被覆 `φ : Y → Z` を通すと最小定義体は
**真に小さくなりうる**（原文 `Proposition 1.7` はまさにその現象を使う）。
★同型のときだけ両向きが取れて等号になる。★★**その差が `log-diff_Y − log-diff_Z ≥ 0`
という原文の非対称の出どころ**である。

## ★逸脱の記録（CLAUDE.md の「逸脱」）

★原文は「`X_ℚ` にしか依らない」を **ℚ-スキームとしての同一性**で言うが、
本ファイルは **`X ≅ X'`（同型）で不変**という形で取る。
★★原文の主張は「`X_ℚ ≅ X'_ℚ` なら同じ」であり、我々の `xF : Spec F ⟶ X` は
ℚ-スキームの点なので、同型による不変性がそのまま原文の主張になる。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField

/-- ★★★★★★**最小定義体は同型で不変**。

原文 (GenEll p.8):
> function log-diffX on X(Q) [cf. Definition 1.5, (iii)] depends only on the scheme

★`minField_comp_le`（任意の射で `≤`）を `e.hom` と `e.inv` に当てるだけ。 -/
theorem minField_iso (F : Type) [Field F] {X X' : Scheme} (e : X ≅ X')
    (xF : Spec (CommRingCat.of F) ⟶ X) :
    minField F (xF ≫ e.hom) = minField F xF := by
  refine le_antisymm (minField_comp_le F e.hom xF) ?_
  have h := minField_comp_le F e.inv (xF ≫ e.hom)
  rwa [Category.assoc, e.hom_inv_id, Category.comp_id] at h

/-- ★★★★★★★**[GenEll] Remark 1.5.1 の (A)** —— `log-diff` は同型で不変。

原文 (GenEll p.8):
> function log-diffX on X(Q) [cf. Definition 1.5, (iii)] depends only on the scheme

★★★**これが「`log-diff_X` は `X_ℚ` にしか依らない」の中身**である。
★spreading out は要らない——`log-diff` は点の最小定義体だけで決まるからである。

★★対比: `log-cond_D` のほうは ℤ-モデルに依り、BD-類でしか一致しない（(B)）。
その計算部分は `LogCondSigma.lean` の `abs_logCond_sub_le_sum_log` にあるが、
★**2 つのモデルが Σ の外で一致すること**（spreading out）は未実装である。 -/
theorem logDiff_minField_iso (F : Type) [Field F] [NumberField F]
    {X X' : Scheme} (e : X ≅ X') (xF : Spec (CommRingCat.of F) ⟶ X) :
    haveI := NumberField.of_subfield (minField F xF)
    haveI := NumberField.of_subfield (minField F (xF ≫ e.hom))
    logDiffOfField (minField F (xF ≫ e.hom)) = logDiffOfField (minField F xF) := by
  haveI := NumberField.of_subfield (minField F xF)
  haveI := NumberField.of_subfield (minField F (xF ≫ e.hom))
  exact logDiffOfField_congr (RingEquiv.subfieldCongr (minField_iso F e xF))

/-! ### ★出典の紐付け(`.src`)

★★**項目全体の `.src` は置かない。** `Remark 1.5.1` は (A)(B) の 2 主張からなり、
(B) は算術的 spreading out（`ResearchPaper/mathlib-gap.json` の
`spreading-out-over-base`）が残っている。本ファイルは **(A) だけ**を閉じている。 -/

def minField_iso.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Remark 1.5.1(最小定義体は同型で不変)",
    sectionId := "genell-rem-1-5-1" }

def logDiff_minField_iso.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Remark 1.5.1, (A)(log-diff は同型で不変——log-cond の側は含まない)",
    sectionId := "genell-rem-1-5-1" }

def logDiff_minField_iso.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "minField_comp_le(任意の射で最小定義体は小さくなる)"
      (.inProject "ABC3" "ABC3.Found.GenEll.minField_comp_le") 8,
    .citation "[ABC3]" "logDiffOfField_congr(同型な数体は同じ log-diff)"
      (.inProject "ABC3" "ABC3.Found.GenEll.logDiffOfField_congr") 8,
    .implicitStep
      ("★★★Remark 1.5.1 は 2 主張である。(A) log-diff は X_ℚ にしか依らない" ++
       "——本ファイル。(B) log-cond の BD-類は (X_ℚ, D_ℚ) にしか依らない" ++
       "——計算部分は LogCondSigma.lean の abs_logCond_sub_le_sum_log にあるが、" ++
       "2 つのモデルが Σ の外で一致すること(算術的 spreading out)は未実装") 8,
    .implicitStep
      ("★★(A) に spreading out が要らないのは log-diff が点の最小定義体だけで" ++
       "決まるからである(Definition 1.5, (iii))。因子も計量も見ないので" ++
       "ℤ-モデルの取り方が入り込む余地がない") 8,
    .implicitStep
      ("★対比: minField_comp_le は不等号であり、被覆を通すと最小定義体は真に" ++
       "小さくなりうる(Proposition 1.7 はその現象を使う)。同型のときだけ両向きが" ++
       "取れて等号になる") 8 ]

end ABC3.Found.GenEll
