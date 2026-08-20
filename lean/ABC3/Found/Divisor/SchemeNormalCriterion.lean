/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.NormFunctor
import Mathlib.AlgebraicGeometry.AffineScheme
import Mathlib.RingTheory.IntegralClosure.IntegrallyClosed

/-!
# 正規性のアフィン局所判定(鎖 `normalize` の `normalization-proper-normal` の前段)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.109。

原文 (FrdI p.109):
> a proper normal variety], then let us write DL for the set of prime divisors of V [L]

## ★★何のためか

`SchemeWeil.lean` の `IsNormalScheme X`(各茎が整閉)を**アフィン局所**で判定する。
★`V[L]`(`NormFunctor.lean`)の正規性は、アフィン開の切断が**整閉包**であること
(`Scheme.Hom.normalizationObjIso`)から出るので、この橋が要る。

## ★★★中身は 2 本

| 段 | mathlib の在庫 |
|---|---|
| アフィン開の茎は切断の局所化 | `IsAffineOpen.isLocalization_stalk` |
| 整閉性は局所化で保たれる | `isIntegrallyClosed_of_isLocalization` |

★`primeCompl ≤ R⁰`(`Ideal.primeCompl_le_nonZeroDivisors`)が整域であることを使う所。

★★**実装上の注意**: 茎への `Algebra` 構造は
`TopCat.Presheaf.algebra_section_stalk` で**明示に入れる**必要がある
(`x : X` の形では instance 探索が届かない)。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `isNormalScheme_of_affineOpens` | ★★**アフィン開の切断が整閉なら `X` は正規** |
-/

namespace ABC3.Found.Divisor

open AlgebraicGeometry CategoryTheory

universe u

/-- ★★★**正規性のアフィン局所判定** —— すべてのアフィン開の切断が整閉なら、
`X` の各茎は整閉。

★茎はアフィン開の切断の局所化(`IsAffineOpen.isLocalization_stalk`)であり、
整閉性は局所化で保たれる(`isIntegrallyClosed_of_isLocalization`)。 -/
theorem isNormalScheme_of_affineOpens (X : Scheme.{u}) [IsIntegral X]
    (h : ∀ U : X.affineOpens, IsIntegrallyClosed Γ(X, U.1)) :
    IsNormalScheme X := by
  intro x
  obtain ⟨_, ⟨U, hU, rfl⟩, hxU, -⟩ :=
    X.isBasis_affineOpens.exists_subset_of_mem_open (Set.mem_univ x) isOpen_univ
  haveI : Nonempty U := ⟨⟨x, hxU⟩⟩
  haveI := h ⟨U, hU⟩
  haveI : IsDomain Γ(X, U) := inferInstance
  letI : Algebra Γ(X, U) (X.presheaf.stalk ((⟨x, hxU⟩ : U) : X)) :=
    TopCat.Presheaf.algebra_section_stalk X.presheaf ⟨x, hxU⟩
  haveI := hU.isLocalization_stalk ⟨x, hxU⟩
  exact isIntegrallyClosed_of_isLocalization (R := Γ(X, U))
    (X.presheaf.stalk ((⟨x, hxU⟩ : U) : X)) (hU.primeIdealOf ⟨x, hxU⟩).asIdeal.primeCompl
    (Ideal.primeCompl_le_nonZeroDivisors _)

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Example 6.1` の「`V[L]` は正規」を出すためのアフィン局所判定。 -/
def isNormalScheme_of_affineOpens.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — 正規性のアフィン局所判定(V[L] が正規であることの前段)",
    sectionId := "frdi-example-6-1" }

end ABC3.Found.Divisor
