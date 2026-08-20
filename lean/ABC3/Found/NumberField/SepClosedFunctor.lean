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
noncomputable def compObj (L : FinSub K₁ Ω) : FinSub K₂ Ω where
  toIF := IntermediateField.adjoin K₂ (L.toIF : Set Ω)
  fin := by
    haveI := L.fin
    exact finiteDimensional_adjoin_coe (K₂ := K₂) L.toIF (primElt L.toIF)
      (primElt_adjoin L.toIF)

@[simp] theorem compObj_toIF (L : FinSub K₁ Ω) :
    (compObj (K₂ := K₂) L).toIF = IntermediateField.adjoin K₂ (L.toIF : Set Ω) := rfl

/-- ★★★**次数は保たれる**。 -/
theorem finrank_compObj (hsc : separableClosure K₁ K₂ = ⊥) (L : FinSub K₁ Ω) :
    Module.finrank K₂ (compObj (K₂ := K₂) L).toIF = Module.finrank K₁ L.toIF := by
  haveI := L.fin
  exact finrank_adjoin_coe_eq_of_separableClosure_eq_bot hsc L.toIF

/-! ## ★2. 射の側 -/

/-- ★★**`f : L ⟶ M` が定める `K₂(L) →ₐ[K₂] Ω`**。 -/
noncomputable def compHomToOmega (hsc : separableClosure K₁ K₂ = ⊥)
    {L M : FinSub K₁ Ω} (f : L ⟶ M) :
    (compObj (K₂ := K₂) L).toIF →ₐ[K₂] Ω :=
  haveI := L.fin
  haveI := M.fin
  (extendAlgHom hsc L.toIF ((M.toIF.val).comp (FinSub.hom f)) (primElt L.toIF)).comp
    (IntermediateField.equivOfEq (adjoin_coe_eq_adjoin_simple (K₂ := K₂) L.toIF
      (primElt L.toIF) (primElt_adjoin L.toIF))).toAlgHom

/-- ★`L` の元の上では `f` そのもの。 -/
theorem compHomToOmega_apply (hsc : separableClosure K₁ K₂ = ⊥)
    {L M : FinSub K₁ Ω} (f : L ⟶ M) (x : Ω) (hx : x ∈ L.toIF) :
    compHomToOmega hsc f ⟨x, IntermediateField.subset_adjoin K₂ _ hx⟩
      = ((FinSub.hom f) ⟨x, hx⟩ : Ω) := by
  haveI := L.fin
  haveI := M.fin
  have heq := adjoin_coe_eq_adjoin_simple (K₂ := K₂) L.toIF (primElt L.toIF)
    (primElt_adjoin L.toIF)
  have h2 : (IntermediateField.equivOfEq heq)
      ⟨x, IntermediateField.subset_adjoin K₂ _ hx⟩
      = ⟨x, le_adjoin_K₂ L.toIF (primElt L.toIF) (primElt_adjoin L.toIF) hx⟩ := rfl
  show extendAlgHom hsc L.toIF ((M.toIF.val).comp (FinSub.hom f)) (primElt L.toIF)
      ((IntermediateField.equivOfEq heq)
        ⟨x, IntermediateField.subset_adjoin K₂ _ hx⟩) = _
  rw [h2]
  exact extendAlgHom_restrict hsc L.toIF ((M.toIF.val).comp (FinSub.hom f))
    (primElt L.toIF) (primElt_adjoin L.toIF) x hx

/-- ★★**抽象版** —— `adjoin F s` からの代数射の像は、生成元の像で決まる。

★★`compHomToOmega` を展開させないために、射を**抽象のまま**扱うのが要点である
(具体形のまま `adjoin_induction` を回すと `whnf` が落ちる)。 -/
theorem algHom_mem_of_mem_adjoin {F E : Type u} [Field F] [Field E] [Algebra F E]
    (s : Set E) (g : (IntermediateField.adjoin F s) →ₐ[F] E) (T : IntermediateField F E)
    (h : ∀ x (hx : x ∈ s), g ⟨x, IntermediateField.subset_adjoin F s hx⟩ ∈ T)
    (x : E) (hx : x ∈ IntermediateField.adjoin F s) : g ⟨x, hx⟩ ∈ T := by
  refine IntermediateField.adjoin_induction F
    (p := fun y hy => g ⟨y, hy⟩ ∈ T) h ?_ ?_ ?_ ?_ hx
  · intro c
    have h1 : g ⟨algebraMap F E c, IntermediateField.algebraMap_mem _ _⟩ = algebraMap F E c :=
      AlgHom.commutes _ c
    rw [h1]
    exact T.algebraMap_mem c
  · intro y z hy hz hy' hz'
    have h1 : g ⟨y + z, add_mem hy hz⟩ = g ⟨y, hy⟩ + g ⟨z, hz⟩ :=
      map_add g ⟨y, hy⟩ ⟨z, hz⟩
    rw [h1]
    exact add_mem hy' hz'
  · intro y hy hy'
    have h1 : g ⟨y⁻¹, inv_mem hy⟩ = (g ⟨y, hy⟩)⁻¹ := map_inv₀ g ⟨y, hy⟩
    rw [h1]
    exact inv_mem hy'
  · intro y z hy hz hy' hz'
    have h1 : g ⟨y * z, mul_mem hy hz⟩ = g ⟨y, hy⟩ * g ⟨z, hz⟩ :=
      map_mul g ⟨y, hy⟩ ⟨z, hz⟩
    rw [h1]
    exact mul_mem hy' hz'

/-- ★★★**像は `K₂(M)` に入る** —— 生成元での値が `M` に入るから。 -/
theorem compHomToOmega_mem (hsc : separableClosure K₁ K₂ = ⊥)
    {L M : FinSub K₁ Ω} (f : L ⟶ M) (x : Ω)
    (hx : x ∈ IntermediateField.adjoin K₂ (L.toIF : Set Ω)) :
    compHomToOmega hsc f ⟨x, hx⟩
      ∈ IntermediateField.adjoin K₂ (M.toIF : Set Ω) := by
  refine algHom_mem_of_mem_adjoin (L.toIF : Set Ω) (compHomToOmega hsc f) _ ?_ x hx
  intro y hy
  show compHomToOmega hsc f ⟨y, IntermediateField.subset_adjoin K₂ _ hy⟩ ∈ _
  rw [compHomToOmega_apply hsc f y hy]
  exact IntermediateField.subset_adjoin K₂ _ ((FinSub.hom f) ⟨y, hy⟩).2

/-- ★★★★**射の側** —— `f : L ⟶ M` が定める `K₂(L) →ₐ[K₂] K₂(M)`。 -/
noncomputable def compMap (hsc : separableClosure K₁ K₂ = ⊥)
    {L M : FinSub K₁ Ω} (f : L ⟶ M) :
    (compObj (K₂ := K₂) L).toIF →ₐ[K₂] (compObj (K₂ := K₂) M).toIF :=
  AlgHom.codRestrict (compHomToOmega hsc f) ((compObj (K₂ := K₂) M).toIF.toSubalgebra)
    (fun x => compHomToOmega_mem hsc f x.1 x.2)

@[simp] theorem compMap_coe (hsc : separableClosure K₁ K₂ = ⊥)
    {L M : FinSub K₁ Ω} (f : L ⟶ M) (x : (compObj (K₂ := K₂) L).toIF) :
    ((compMap hsc f x : (compObj (K₂ := K₂) M).toIF) : Ω) = compHomToOmega hsc f x := rfl

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Theorem 6.2, (i)` の `L ↦ K₂(L)`。 -/
def compObj.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — L ↦ K₂(L) は次数を保つ",
    sectionId := "frdi-thm-6-2" }

end ABC3.Found.NF
