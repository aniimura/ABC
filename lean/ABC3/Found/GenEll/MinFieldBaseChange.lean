/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.MinFieldCovering

/-!
# [GenEll] Definition 1.5, (i) —— **最小定義体は底変換で対応する**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.8。

原文 (GenEll p.8):
> one obtains a well-deﬁned minimal ﬁeld of deﬁnition Fmin ⊆ F of x.

## ★★何のために要るか

`Definition 1.5, (iii)(iv)` は `log-diff` / `log-cond` を
**最小定義体で測る**ことで定義する（`LogDiffPoint.lean`）。
★それが `X(ℚ̄)` の関数として well-defined であるには、

> **同じ点を別の定義体で表しても、最小定義体が対応する**

ことが要る。本ファイルがそれを取る。

## ★★★中身は「像点が動かない」ことだけ

`MinField.lean` は `F_min` を **mathlib の `SpecToEquivOfField`**

    (Spec K ⟶ X) ≃ Σ ξ : X, (κ(ξ) ⟶ K)

の像として作っている。★底変換 `Spec K → Spec F → X` を этой 全単射で読むと

    ⟨ξ, κ(ξ) → F⟩  ↦  ⟨ξ, κ(ξ) → F → K⟩

であり、**像点 `ξ` は変わらない**（`minFieldData_baseChange`）。
★★したがって `K_min` は `F_min` の `F → K` による像である（`minField_baseChange`）。

★★★**極小性の議論も scheme 論的像の理論も要らなかった** ——
全単射の単射性を 1 回使うだけである。`MinField.lean` が
「またしても正面から作る必要が無かった」と書いたのと同じ機構である。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory

/-- ★★★★**底変換しても像点は変わらない** —— 最小定義体の対応の中身。

原文 (GenEll p.8):
> one obtains a well-deﬁned minimal ﬁeld of deﬁnition Fmin ⊆ F of x.

★`SpecToEquivOfField` の**単射性を 1 回**使うだけである。 -/
theorem minFieldData_baseChange (F K : Type) [Field F] [Field K] [Algebra F K] {X : Scheme}
    (xF : Spec (CommRingCat.of F) ⟶ X) :
    minFieldData K (Spec.map (CommRingCat.ofHom (algebraMap F K)) ≫ xF)
      = ⟨(minFieldData F xF).1,
         (minFieldData F xF).2 ≫ CommRingCat.ofHom (algebraMap F K)⟩ := by
  apply (Scheme.SpecToEquivOfField K X).symm.injective
  rw [show (Scheme.SpecToEquivOfField K X).symm
        (minFieldData K (Spec.map (CommRingCat.ofHom (algebraMap F K)) ≫ xF))
      = Spec.map (CommRingCat.ofHom (algebraMap F K)) ≫ xF from
    (Scheme.SpecToEquivOfField K X).symm_apply_apply _]
  show _ = Spec.map ((minFieldData F xF).2 ≫ CommRingCat.ofHom (algebraMap F K))
      ≫ X.fromSpecResidueField (minFieldData F xF).1
  rw [Spec.map_comp, Category.assoc, specMap_minFieldData_comp]

/-- ★★★★★**最小定義体は底変換で対応する** —— `K_min` は `F_min` の像である。

原文 (GenEll p.8):
> one obtains a well-deﬁned minimal ﬁeld of deﬁnition Fmin ⊆ F of x.

★★★これが `Definition 1.5, (iii)(iv)` の
「`log-diff` / `log-cond` を最小定義体で測る」が `X(ℚ̄)` 上で well-defined である
ことの**幾何の側**である。★残るのは代数の側
（同型な数体は同じ差積次数をもつ）だけになった。 -/
theorem minField_baseChange (F K : Type) [Field F] [Field K] [Algebra F K] {X : Scheme}
    (xF : Spec (CommRingCat.of F) ⟶ X) :
    minField K (Spec.map (CommRingCat.ofHom (algebraMap F K)) ≫ xF)
      = (minField F xF).map (algebraMap F K) := by
  rw [minField, minFieldData_baseChange F K xF]
  ext y
  constructor
  · rintro ⟨a, rfl⟩
    exact ⟨(minFieldData F xF).2.hom a, ⟨a, rfl⟩, rfl⟩
  · rintro ⟨_, ⟨a, rfl⟩, rfl⟩
    exact ⟨a, rfl⟩

/-- ★★**最小定義体であることは底変換で保たれない** ——
`F` が最小でも `K ⊋ F` では最小でない。

★これは `minField_baseChange` の直接の帰結である:
`K_min = F_min` の像 ⊆ `F` の像 ⊊ `K`。
★★**`log-diff` が底変換で増える**（`logDiffAt_le_baseChange`）ことと表裏である。 -/
theorem not_isMinimalFieldOfDefinition_baseChange (F K : Type) [Field F] [Field K]
    [Algebra F K] {X : Scheme} (xF : Spec (CommRingCat.of F) ⟶ X)
    (h : ((minField F xF).map (algebraMap F K)) ≠ ⊤) :
    ¬ IsMinimalFieldOfDefinition K (Spec.map (CommRingCat.ofHom (algebraMap F K)) ≫ xF) := by
  intro hmin
  exact h (by rw [← minField_baseChange F K xF]; exact hmin)

/-! ### ★出典の紐付け(`.src`) -/

def minFieldData_baseChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (i)(底変換しても像点は変わらない)",
    sectionId := "genell-def-1-5" }

def minField_baseChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (i)(最小定義体は底変換で対応する)",
    sectionId := "genell-def-1-5" }

def minField_baseChange.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Scheme.SpecToEquivOfField((Spec K ⟶ X) ≃ Σ ξ, (κ(ξ) ⟶ K))"
      (.inMathlib "AlgebraicGeometry.Scheme.SpecToEquivOfField") 8,
    .citation "[ABC3]" "specMap_minFieldData_comp(F_min が実際に定義体であること)"
      (.inProject "ABC3" "ABC3.Found.GenEll.specMap_minFieldData_comp") 8,
    .implicitStep
      ("★★これで Definition 1.5, (iii)(iv) の well-defined 性の**幾何の側**は取れた。" ++
       "☆残るのは代数の側 —— 同型な数体は同じ差積次数をもつこと") 8 ]

def not_isMinimalFieldOfDefinition_baseChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (i)(最小性は底変換で保たれない)",
    sectionId := "genell-def-1-5" }

end ABC3.Found.GenEll
