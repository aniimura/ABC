/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.NumberField.SepClosedIn
import ABC3.Found.FrdI.Sec6GaloisCat

/-!
# `𝒟₁ → 𝒟₂` の関手(鎖 `sec6items` の `thm62-i-Dfun`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.110。

原文 (FrdI p.110):
> for the group of rational functions on Vi[Li] whose zeroes and poles are supported

## ★★`L ↦ K₂(L)` が関手になる

`K₁` が `K₂` の中で分離閉なら、`L ↦ K₂(L)` は
`FinSub K₁ Ω ⥤ FinSub K₂ Ω` を定める。

| 段 | 根拠(`SepClosedIn.lean`) |
|---|---|
| 対象が有限次元 | `finiteDimensional_adjoin_coe` |
| 次数の保存 | `finrank_adjoin_coe_eq_of_separableClosure_eq_bot` |
| 射 | `extendAlgHom` ＋ `extendAlgHom_restrict` |
| 関手則 | `algHom_adjoin_ext`(延長の一意性) |

## ★本ファイルで閉じること

| 定義 | 中身 |
|---|---|
| `primElt` | 原始元の選択 |
| `compObj` | ★`L ↦ K₂(L)`(対象の側) |
| `finrank_compObj` | ★★**次数は保たれる** |
-/

namespace ABC3.Found.NF

open Polynomial CategoryTheory ABC3.Found.FrdI IntermediateField

universe u

variable {K₁ K₂ Ω : Type u} [Field K₁] [Field K₂] [Field Ω]
  [Algebra K₁ K₂] [Algebra K₁ Ω] [Algebra K₂ Ω] [IsScalarTower K₁ K₂ Ω]
  [Algebra.IsSeparable K₁ Ω]

/-- ★原始元の選択。 -/
noncomputable def primElt (L : IntermediateField K₁ Ω) [FiniteDimensional K₁ L] : L :=
  (Field.exists_primitive_element K₁ L).choose

theorem primElt_spec (L : IntermediateField K₁ Ω) [FiniteDimensional K₁ L] :
    K₁⟮(primElt L : L)⟯ = ⊤ :=
  (Field.exists_primitive_element K₁ L).choose_spec

/-- ★原始元は `Ω` の側でも `L` を生成する。 -/
theorem primElt_adjoin (L : IntermediateField K₁ Ω) [FiniteDimensional K₁ L] :
    K₁⟮((primElt L : L) : Ω)⟯ = L := by
  have h := congrArg IntermediateField.lift (primElt_spec L)
  rwa [IntermediateField.lift_adjoin_simple, IntermediateField.lift_top] at h

/-- ★★**対象の側** —— `L ↦ K₂(L)`。 -/
noncomputable def compObj (L : FinSub K₁ Ω) : FinSub K₂ Ω :=
  haveI := L.fin
  { toIF := IntermediateField.adjoin K₂ (L.toIF : Set Ω)
    fin := finiteDimensional_adjoin_coe (K₂ := K₂) L.toIF (primElt L.toIF)
      (primElt_adjoin L.toIF) }

/-- ★★★**次数は保たれる**。 -/
theorem finrank_compObj (hsc : separableClosure K₁ K₂ = ⊥) (L : FinSub K₁ Ω) :
    Module.finrank K₂ (compObj (K₂ := K₂) L).toIF = Module.finrank K₁ L.toIF := by
  haveI := L.fin
  exact finrank_adjoin_coe_eq_of_separableClosure_eq_bot hsc L.toIF

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Theorem 6.2, (i)` の `L ↦ K₂(L)`。 -/
def compObj.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — L ↦ K₂(L) は次数を保つ",
    sectionId := "frdi-thm-6-2" }

end ABC3.Found.NF
