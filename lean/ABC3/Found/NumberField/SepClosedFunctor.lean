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

/-! ## ★3. 関手則 -/

/-- ★★**合成則**(`Ω` の側で) —— 延長の一意性から。 -/
theorem compHomToOmega_comp (hsc : separableClosure K₁ K₂ = ⊥)
    {L M N : FinSub K₁ Ω} (f : L ⟶ M) (g : M ⟶ N) :
    compHomToOmega hsc (f ≫ g) = (compHomToOmega hsc g).comp (compMap hsc f) := by
  refine algHom_adjoin_ext (K₂ := K₂) (L.toIF : Set Ω) ?_
  intro y hy
  show compHomToOmega hsc (f ≫ g) ⟨y, IntermediateField.subset_adjoin K₂ _ hy⟩ = _
  rw [compHomToOmega_apply hsc (f ≫ g) y hy]
  show ((FinSub.hom (f ≫ g)) ⟨y, hy⟩ : Ω)
    = compHomToOmega hsc g (compMap hsc f ⟨y, IntermediateField.subset_adjoin K₂ _ hy⟩)
  have h1 : (compMap hsc f ⟨y, IntermediateField.subset_adjoin K₂ _ hy⟩
      : (compObj (K₂ := K₂) M).toIF)
      = ⟨((FinSub.hom f) ⟨y, hy⟩ : Ω),
        IntermediateField.subset_adjoin K₂ _ ((FinSub.hom f) ⟨y, hy⟩).2⟩ := by
    refine Subtype.ext ?_
    rw [compMap_coe]
    exact compHomToOmega_apply hsc f y hy
  rw [h1, compHomToOmega_apply hsc g _ ((FinSub.hom f) ⟨y, hy⟩).2]
  rfl

/-- ★★**恒等則**(`Ω` の側で)。 -/
theorem compHomToOmega_id (hsc : separableClosure K₁ K₂ = ⊥) (L : FinSub K₁ Ω) :
    compHomToOmega hsc (𝟙 L) = (compObj (K₂ := K₂) L).toIF.val := by
  refine algHom_adjoin_ext (K₂ := K₂) (L.toIF : Set Ω) ?_
  intro y hy
  show compHomToOmega hsc (𝟙 L) ⟨y, IntermediateField.subset_adjoin K₂ _ hy⟩ = _
  rw [compHomToOmega_apply hsc (𝟙 L) y hy]
  rfl

/-- ★★★★★**`L ↦ K₂(L)` は関手** —— `𝒟₁ → 𝒟₂`(原文 `Theorem 6.2, (i)`)。 -/
noncomputable def compFunctor (hsc : separableClosure K₁ K₂ = ⊥) :
    FinSub K₁ Ω ⥤ FinSub K₂ Ω where
  obj := compObj
  map f := compMap hsc f
  map_id L := by
    refine AlgHom.ext fun x => Subtype.ext ?_
    rw [compMap_coe]
    show compHomToOmega hsc (𝟙 L) x = _
    rw [compHomToOmega_id hsc L]
    rfl
  map_comp f g := by
    refine AlgHom.ext fun x => Subtype.ext ?_
    rw [compMap_coe]
    show compHomToOmega hsc (f ≫ g) x = _
    rw [compHomToOmega_comp hsc f g]
    rfl

/-! ## ★4. 埋め込みに沿った移送(`Theorem 6.2, (i)` の仮定 (b))

★原文の仮定 (b)「`K₁ → K₂ → K̄₂` が `K̄₁` を経由する」は、
**`ι : K̄₁ →ₐ[K₁] K̄₂` が取れる**ということである。
★このとき `𝒟₁ → 𝒟₂` は「`ι` に沿った移送」と `compFunctor` の**合成**に分かれる。 -/

section Transport

variable {K Kbar₁ Kbar₂ : Type u} [Field K] [Field Kbar₁] [Field Kbar₂]
  [Algebra K Kbar₁] [Algebra K Kbar₂] (ι : Kbar₁ →ₐ[K] Kbar₂)

/-- ★**対象の移送** —— `L ↦ ι(L)`。 -/
noncomputable def mapObj (L : FinSub K Kbar₁) : FinSub K Kbar₂ where
  toIF := L.toIF.map ι
  fin := by
    haveI := L.fin
    exact LinearEquiv.finiteDimensional (IntermediateField.equivMap L.toIF ι).toLinearEquiv

@[simp] theorem mapObj_toIF (L : FinSub K Kbar₁) :
    (mapObj ι L).toIF = L.toIF.map ι := rfl

/-- ★**射の移送** —— 同型 `L ≃ ι(L)` で共役する。 -/
noncomputable def mapHomF {L M : FinSub K Kbar₁} (f : L ⟶ M) :
    (mapObj ι L).toIF →ₐ[K] (mapObj ι M).toIF :=
  ((IntermediateField.equivMap M.toIF ι).toAlgHom).comp
    ((FinSub.hom f).comp (IntermediateField.equivMap L.toIF ι).symm.toAlgHom)

/-- ★★★**`L ↦ ι(L)` は関手** —— `FinSub K K̄₁ ⥤ FinSub K K̄₂`。 -/
noncomputable def mapFinSub : FinSub K Kbar₁ ⥤ FinSub K Kbar₂ where
  obj := mapObj ι
  map f := mapHomF ι f
  map_id L := by
    refine AlgHom.ext fun x => ?_
    show (IntermediateField.equivMap L.toIF ι)
      ((FinSub.hom (𝟙 L)) ((IntermediateField.equivMap L.toIF ι).symm x)) = x
    rw [FinSub.hom_id]
    show (IntermediateField.equivMap L.toIF ι)
      ((IntermediateField.equivMap L.toIF ι).symm x) = x
    rw [AlgEquiv.apply_symm_apply]
  map_comp f g := by
    refine AlgHom.ext fun x => ?_
    show (IntermediateField.equivMap _ ι)
      ((FinSub.hom (f ≫ g)) ((IntermediateField.equivMap _ ι).symm x)) = _
    rw [FinSub.hom_comp]
    show (IntermediateField.equivMap _ ι)
      ((FinSub.hom g) ((FinSub.hom f) ((IntermediateField.equivMap _ ι).symm x))) = _
    show _ = (IntermediateField.equivMap _ ι)
      ((FinSub.hom g) ((IntermediateField.equivMap _ ι).symm
        ((IntermediateField.equivMap _ ι) ((FinSub.hom f)
          ((IntermediateField.equivMap _ ι).symm x)))))
    rw [AlgEquiv.symm_apply_apply]

end Transport

/-! ## ★5. 原文の形の `𝒟₁ → 𝒟₂` -/

/-- ★★★★★★**原文 `Theorem 6.2, (i)` の関手 `𝒟₁ → 𝒟₂`**。

★仮定 (b)(`K̄₁ → K̄₂` の埋め込み `ι`)と (c)(`K₁` が `K₂` の中で分離閉)から、
`𝒟ᵢ = B(Gᵢ)⁰ = (FinSub Kᵢ K̄ᵢ)ᵒᵖ` の間の関手が
**「`ι` に沿った移送」と「`L ↦ K₂(L)`」の合成の反対**として出る。 -/
noncomputable def dFunctor {K₁ K₂ Kbar₁ Kbar₂ : Type u} [Field K₁] [Field K₂]
    [Field Kbar₁] [Field Kbar₂] [Algebra K₁ K₂] [Algebra K₁ Kbar₁] [Algebra K₁ Kbar₂]
    [Algebra K₂ Kbar₂] [IsScalarTower K₁ K₂ Kbar₂] [Algebra.IsSeparable K₁ Kbar₂]
    (ι : Kbar₁ →ₐ[K₁] Kbar₂) (hsc : separableClosure K₁ K₂ = ⊥) :
    (FinSub K₁ Kbar₁)ᵒᵖ ⥤ (FinSub K₂ Kbar₂)ᵒᵖ :=
  (mapFinSub ι ⋙ compFunctor hsc).op

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Theorem 6.2, (i)` の `L ↦ K₂(L)`。 -/
def compObj.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — L ↦ K₂(L) は次数を保つ",
    sectionId := "frdi-thm-6-2" }

/-- ★★★★★locator —— `Theorem 6.2, (i)` の関手 `𝒟₁ → 𝒟₂`。 -/
def compFunctor.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — ψ が 𝒟₁ → 𝒟₂ を定める",
    sectionId := "frdi-thm-6-2" }

/-- ★★★★★★locator —— `Theorem 6.2, (i)` の仮定 (b)(c) から出る `𝒟₁ → 𝒟₂`。 -/
def dFunctor.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — 体の包含 K₁ → K₂ が誘導する 𝒟₁ → 𝒟₂",
    sectionId := "frdi-thm-6-2" }

end ABC3.Found.NF
