/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm51Span

/-!
# [FrdI] Theorem 5.1 —— span の類は well-defined で、Frobenius-trivial なら 0

`Thm51Span.lean` で「span 1 個の判定条件」を作った。本ファイルではその上に

1. **類 `spanCls` は span の取り方に依らない**(`Φ^birat` を法として)
   —— これが原文の「写像 `Φ(A) × Φ(A) → Pic_𝒞(A)` は `Pic_Φ(A)` を経由する」
2. **Frobenius-trivial な対象では類が `2` 倍不変、ゆえに `0`**
   —— これが原文の「`d · ξ = ξ` for all `d`、`d = 2` で `ξ = 0`」

を積む。★2 の帰結が `Theorem 5.1, (iii)` の中心である:
**base-isomorphic な Frobenius-trivial 対象は同型で、すべて Aut-ample**。

原文 (FrdI p.99):
> morphisms of Frobenius type of arbitrary prescribed Frobenius degree] corresponds

原文 (FrdI p.99):
> and that all Frobenius-trivial objects of C are Aut-ample. In light of these
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-! ## ★1. 類の well-defined 性 -/

/-- ★span を `Z` へ引き戻したときの差の membership。 -/
theorem span_pullback_mem (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {Z X A' : C} (u : Z ⟶ X) (hu : IsPreStep P u) (φ : X ⟶ A') (hsφ : IsPreStep P φ)
    {X' : C} (u' : Z ⟶ X') (hu' : IsPreStep P u') (φ' : X' ⟶ A') (hsφ' : IsPreStep P φ')
    (hb : P.Base (u ≫ φ) = P.Base (u' ≫ φ')) :
    toGp _ (P.Div (u' ≫ φ')) - toGp _ (P.Div (u ≫ φ)) ∈ phiBiratAt P G Z := by
  have hs1 : IsPreStep P (u ≫ φ) := IsPreStep.comp P hu hsφ
  have hs2 : IsPreStep P (u' ≫ φ') := IsPreStep.comp P hu' hsφ'
  have hmem := mem_phiBiratAt_of_preStepPair G (u ≫ φ) (u' ≫ φ')
    (prop_1_4_i P _ (fun Y _ => hiso Y)) hs1
    (prop_1_4_i P _ (fun Y _ => hiso Y)) hs2 hb
  rw [← spanCls_eq_sliceDivGpOf _ _ _ hs2.1] at hmem
  have htr := mem_phiBiratAt_transport G (u ≫ φ) (prop_1_4_i P _ (fun Y _ => hiso Y)) hs1
    (spanCls (u ≫ φ) hs1.2 (u' ≫ φ'))
  rw [gpMap_base_spanCls] at htr
  exact htr.mpr hmem

/-- ★span を pre-step `u` で引き戻したときの類の押し出し(同じ頂点側)。 -/
theorem gpMap_base_comp_spanCls {Z X A' A : C} (u : Z ⟶ X) (φ : X ⟶ A') (ψ : X ⟶ A)
    (hsφ : IsPreStep P φ) (hsψ : IsPreStep P ψ) :
    gpMap _ (Φ.map (P.Base (u ≫ φ))) (spanCls φ hsφ.2 ψ)
      = toGp _ (P.Div (u ≫ ψ)) - toGp _ (P.Div (u ≫ φ)) := by
  rw [P.Base_comp, ← gpMap_phi_comp, gpMap_base_spanCls,
    Div_comp_preStep _ _ hsψ.1, Div_comp_preStep _ _ hsφ.1,
    toGp_add, toGp_add, ← gpMap_toGp, ← gpMap_toGp, map_sub]
  abel

/-- ★同じく、もう一方の頂点側。 -/
theorem gpMap_base_comp_spanCls' {Z X X' A' A : C} (u : Z ⟶ X) (u' : Z ⟶ X')
    (φ : X ⟶ A') (φ' : X' ⟶ A') (ψ' : X' ⟶ A)
    (hsφ' : IsPreStep P φ') (hsψ' : IsPreStep P ψ')
    (hu : IsPreStep P u)
    (hspan : P.Base φ ≫ @inv _ _ _ _ (P.Base φ') hsφ'.2
      = @inv _ _ _ _ (P.Base u) hu.2 ≫ P.Base u') :
    gpMap _ (Φ.map (P.Base (u ≫ φ))) (spanCls φ' hsφ'.2 ψ')
      = toGp _ (P.Div (u' ≫ ψ')) - toGp _ (P.Div (u' ≫ φ')) := by
  haveI hbφ' : IsIso (P.Base φ') := hsφ'.2
  haveI hbu : IsIso (P.Base u) := hu.2
  have hinner : gpMap _ (Φ.map (P.Base φ)) (spanCls φ' hsφ'.2 ψ')
      = gpMap _ (Φ.map (inv (P.Base u)))
        (gpMap _ (Φ.map (P.Base u')) (toGp _ (P.Div ψ') - toGp _ (P.Div φ'))) := by
    rw [spanCls_eq, gpMap_phi_comp, hspan, ← gpMap_phi_comp]
  rw [P.Base_comp, ← gpMap_phi_comp, hinner, gpMap_phi_inv_left,
    Div_comp_preStep _ _ hsψ'.1, Div_comp_preStep _ _ hsφ'.1,
    toGp_add, toGp_add, ← gpMap_toGp, ← gpMap_toGp, map_sub]
  abel

/-- ★★★★**類は span の取り方に依らない**(`Φ^birat` を法として)。

原文 (FrdI p.97):
> Now I claim that the map Φ(A) × Φ(A) → PicC(A) of the above discussion factors
-/
theorem span_class_sub_mem (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {X X' A' A : C} (φ : X ⟶ A') (ψ : X ⟶ A) (φ' : X' ⟶ A') (ψ' : X' ⟶ A)
    (hsφ : IsPreStep P φ) (hsψ : IsPreStep P ψ)
    (hsφ' : IsPreStep P φ') (hsψ' : IsPreStep P ψ')
    (hα : @inv _ _ _ _ (P.Base φ) hsφ.2 ≫ P.Base ψ
        = @inv _ _ _ _ (P.Base φ') hsφ'.2 ≫ P.Base ψ') :
    spanCls φ hsφ.2 ψ - spanCls φ' hsφ'.2 ψ' ∈ phiBiratAt P G A' := by
  haveI hbφ : IsIso (P.Base φ) := hsφ.2
  haveI hbψ : IsIso (P.Base ψ) := hsψ.2
  haveI hbφ' : IsIso (P.Base φ') := hsφ'.2
  haveI hbψ' : IsIso (P.Base ψ') := hsψ'.2
  obtain ⟨Z, u, u', hu, hu', hspan⟩ :=
    G.core.preStepSpan X X' (P.Base φ ≫ inv (P.Base φ')) inferInstance
  haveI hbu : IsIso (P.Base u) := hu.2
  haveI hbu' : IsIso (P.Base u') := hu'.2
  have hkey : P.Base u ≫ P.Base φ ≫ inv (P.Base φ') = P.Base u' := by
    rw [hspan, ← Category.assoc, IsIso.hom_inv_id, Category.id_comp]
  have hbφz : P.Base (u ≫ φ) = P.Base (u' ≫ φ') := by
    rw [P.Base_comp, P.Base_comp, ← hkey, Category.assoc, Category.assoc,
      IsIso.inv_hom_id, Category.comp_id]
  have hbψz : P.Base (u ≫ ψ) = P.Base (u' ≫ ψ') := by
    rw [P.Base_comp, P.Base_comp, ← hkey, Category.assoc, Category.assoc, ← hα,
      ← Category.assoc (P.Base φ), IsIso.hom_inv_id, Category.id_comp]
  have hE₁ := span_pullback_mem G hiso u hu φ hsφ u' hu' φ' hsφ' hbφz
  have hE₂ := span_pullback_mem G hiso u hu ψ hsψ u' hu' ψ' hsψ' hbψz
  have hs1 : IsPreStep P (u ≫ φ) := IsPreStep.comp P hu hsφ
  have htr := mem_phiBiratAt_transport G (u ≫ φ) (prop_1_4_i P _ (fun Y _ => hiso Y)) hs1
    (spanCls φ hsφ.2 ψ - spanCls φ' hsφ'.2 ψ')
  refine htr.mp ?_
  rw [map_sub, gpMap_base_comp_spanCls u φ ψ hsφ hsψ,
    gpMap_base_comp_spanCls' u u' φ φ' ψ' hsφ' hsψ' hu hspan]
  have hrw : (toGp (Φ.val (P.toElem.obj Z).base) (P.Div (u ≫ ψ))
        - toGp _ (P.Div (u ≫ φ)))
      - (toGp _ (P.Div (u' ≫ ψ')) - toGp _ (P.Div (u' ≫ φ')))
      = (toGp _ (P.Div (u' ≫ φ')) - toGp _ (P.Div (u ≫ φ)))
        - (toGp _ (P.Div (u' ≫ ψ')) - toGp _ (P.Div (u ≫ ψ))) := by abel
  rw [hrw]
  exact AddSubgroup.sub_mem _ hE₁ hE₂

/-! ## ★2. Frobenius 型の自己射で span をずらす -/

/-- ★★★**Frobenius 型で span をずらす** ——
`φ : X → A'` を pre-step、`κ` を `A'` の base-identity な Frobenius 型自己射、
`μ : X → Y` を同じ次数の Frobenius 型射とすると、
`μ ≫ φ₂ = φ ≫ κ` となる pre-step `φ₂ : Y → A'` が取れる。

★`Definition 1.3, (iv)` の分解と `(ii)` の一意性から。 -/
theorem frobShift (G : Frobenioid P) {X A' A'' Y : C} (φ : X ⟶ A') (hsφ : IsPreStep P φ)
    (κ : A' ⟶ A'') (hκb : IsIso (P.Base κ))
    (μ : X ⟶ Y) (hμ : IsFrobeniusType P μ) (hdeg : P.degFr μ = P.degFr κ) :
    ∃ φ₂ : Y ⟶ A'', IsPreStep P φ₂ ∧ μ ≫ φ₂ = φ ≫ κ := by
  haveI := hκb
  obtain ⟨U, V, γ, β, α, hfac, hγ, hβ, hα⟩ := G.core.arbFactor (φ ≫ κ)
  have hdegα : P.degFr α = 1 := (G.core.pullBackLB α hα).2
  have hdegc : P.degFr (φ ≫ κ) = P.degFr κ := by
    rw [P.degFr_comp, hsφ.1, mul_one]
  have hdegγ : P.degFr γ = P.degFr μ := by
    have h := congrArg P.degFr hfac
    rw [hdegc, P.degFr_comp, P.degFr_comp, hβ.1, hdegα, one_mul, one_mul] at h
    rw [hdeg, h]
  obtain ⟨ρ, hρ, hργ⟩ := G.core.frobDegUniq X U Y γ μ hγ hμ hdegγ
  haveI := hρ
  refine ⟨inv ρ ≫ β ≫ α, ?_, ?_⟩
  · haveI hbφ : IsIso (P.Base φ) := hsφ.2
    haveI hbφκ : IsIso (P.Base (φ ≫ κ)) := by rw [P.Base_comp]; infer_instance
    haveI hbγ : IsIso (P.Base γ) := hγ.2
    haveI hbβ : IsIso (P.Base β) := hβ.2
    haveI hbα : IsIso (P.Base α) := by
      have h := congrArg P.Base hfac
      rw [P.Base_comp, P.Base_comp, P.Base_comp] at h
      haveI h2 : IsIso (P.Base γ ≫ P.Base β ≫ P.Base α) := by
        rw [← h, ← P.Base_comp]; exact hbφκ
      haveI h3 : IsIso (P.Base β ≫ P.Base α) :=
        IsIso.of_isIso_comp_left (P.Base γ) (P.Base β ≫ P.Base α)
      exact IsIso.of_isIso_comp_left (P.Base β) (P.Base α)
    refine IsPreStep.comp P ⟨degFr_of_isIso P (inv ρ), base_isIso_of_iso (asIso ρ).symm⟩ ?_
    exact IsPreStep.comp P hβ ⟨hdegα, hbα⟩
  · rw [← hργ, Category.assoc, ← Category.assoc ρ, IsIso.hom_inv_id, Category.id_comp]
    exact hfac.symm

/-- ★ずらした pre-step の底と零因子。 -/
theorem frobShift_data {X A' A'' Y : C} (φ : X ⟶ A')
    (κ : A' ⟶ A'') (hκd : P.Div κ = 0)
    (μ : X ⟶ Y) (hμd : P.Div μ = 0)
    (φ₂ : Y ⟶ A'') (heq : μ ≫ φ₂ = φ ≫ κ) :
    P.Base μ ≫ P.Base φ₂ = P.Base φ ≫ P.Base κ ∧
      Φ.map (P.Base μ) (P.Div φ₂) = ((P.degFr κ : ℕ+) : ℕ) • P.Div φ := by
  constructor
  · have h := congrArg P.Base heq
    rw [P.Base_comp, P.Base_comp] at h
    exact h
  · have h := congrArg P.Div heq
    rw [P.Div_comp, P.Div_comp, hκd, hμd, map_zero, smul_zero, add_zero, zero_add] at h
    exact h

/-- ★★★**ずらした span の類は `n` 倍、`α` は同じ**。 -/
theorem spanCls_frobShift {X A' A Y : C} (φ : X ⟶ A') (ψ : X ⟶ A)
    (hsφ : IsPreStep P φ)
    (μ : X ⟶ Y) (hμb : IsIso (P.Base μ))
    (φ₂ : Y ⟶ A') (ψ₂ : Y ⟶ A) (hsφ₂ : IsPreStep P φ₂)
    (n : ℕ)
    (hbφ₂ : P.Base μ ≫ P.Base φ₂ = P.Base φ)
    (hbψ₂ : P.Base μ ≫ P.Base ψ₂ = P.Base ψ)
    (hdφ₂ : Φ.map (P.Base μ) (P.Div φ₂) = n • P.Div φ)
    (hdψ₂ : Φ.map (P.Base μ) (P.Div ψ₂) = n • P.Div ψ) :
    spanCls φ₂ hsφ₂.2 ψ₂ = n • spanCls φ hsφ.2 ψ ∧
      (@inv _ _ _ _ (P.Base φ₂) hsφ₂.2 ≫ P.Base ψ₂
        = @inv _ _ _ _ (P.Base φ) hsφ.2 ≫ P.Base ψ) := by
  haveI := hμb
  haveI hbφ : IsIso (P.Base φ) := hsφ.2
  haveI hbφ2 : IsIso (P.Base φ₂) := hsφ₂.2
  have hbφ₂' : P.Base φ₂ = inv (P.Base μ) ≫ P.Base φ := by
    rw [← hbφ₂, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
  have hinv : inv (P.Base φ₂) = inv (P.Base φ) ≫ P.Base μ := by
    refine IsIso.inv_eq_of_hom_inv_id ?_
    rw [hbφ₂', Category.assoc, ← Category.assoc (P.Base φ), IsIso.hom_inv_id,
      Category.id_comp, IsIso.inv_hom_id]
  constructor
  · rw [spanCls_eq, spanCls_eq, hinv, ← gpMap_phi_comp]
    rw [map_sub, gpMap_toGp, gpMap_toGp, hdφ₂, hdψ₂, toGp_nsmul, toGp_nsmul, ← smul_sub]
    rw [← map_nsmul]
  · rw [hinv, Category.assoc, hbψ₂]

/-! ## ★3. ★★★★`Theorem 5.1, (iii)` の中心 -/

/-- ★★★★**[FrdI] Theorem 5.1, (iii)** ——
Frobenius-trivial な 2 つの対象の底の同型は、**対象の同型に持ち上がる**。

原文 (FrdI p.99):
> that there exists an isomorphism α : A →∼ A such that αD = Base(α). In
-/
theorem frobTrivial_iso_of_baseIso (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {A A' : C} (hA : IsFrobeniusTrivial P A) (hA' : IsFrobeniusTrivial P A')
    (α : (P.toElem.obj A').base ⟶ (P.toElem.obj A).base) [IsIso α] :
    ∃ θ : A' ≅ A, P.Base θ.hom = α := by
  obtain ⟨X, φ, ψ, hsφ, hsψ, hα⟩ := G.core.preStepSpan A' A α inferInstance
  obtain ⟨ζ', hdeg', hprop'⟩ := hA'
  obtain ⟨ζ, hdeg, hprop⟩ := hA
  set κ' : A' ⟶ A' := ((ζ' 2 : End A') : A' ⟶ A') with hκ'def
  set κ : A ⟶ A := ((ζ 2 : End A) : A ⟶ A) with hκdef
  have hκ'b : P.Base κ' = P.Base (𝟙 A') := (hprop' 2).1
  have hκb : P.Base κ = P.Base (𝟙 A) := (hprop 2).1
  have hκ'd : P.Div κ' = 0 := (hprop' 2).2.1.2
  have hκd : P.Div κ = 0 := (hprop 2).2.1.2
  have hκ'deg : P.degFr κ' = 2 := hdeg' 2
  have hκdeg : P.degFr κ = 2 := hdeg 2
  obtain ⟨Y, μ, hμ, hμdeg⟩ := G.core.frobDegSurj X 2
  have hμd : P.Div μ = 0 := hμ.1.2
  haveI hμb : IsIso (P.Base μ) := hμ.2
  obtain ⟨φ₂, hsφ₂, heqφ⟩ := frobShift G φ hsφ κ' hκ'b μ hμ (by rw [hμdeg, hκ'deg])
  obtain ⟨ψ₂, hsψ₂, heqψ⟩ := frobShift G ψ hsψ κ hκb μ hμ (by rw [hμdeg, hκdeg])
  obtain ⟨hbφ₂, hdφ₂⟩ := frobShift_data φ κ' hκ'b hκ'd μ hμd φ₂ heqφ
  obtain ⟨hbψ₂, hdψ₂⟩ := frobShift_data ψ κ hκb hκd μ hμd ψ₂ heqψ
  rw [hκ'deg] at hdφ₂
  rw [hκdeg] at hdψ₂
  obtain ⟨hcls, hαeq⟩ := spanCls_frobShift φ ψ hsφ μ hμb φ₂ ψ₂ hsφ₂ 2
    hbφ₂ hbψ₂ hdφ₂ hdψ₂
  have hsub := span_class_sub_mem G hiso φ ψ φ₂ ψ₂ hsφ hsψ hsφ₂ hsψ₂ hαeq.symm
  rw [hcls] at hsub
  have hc : spanCls φ hsφ.2 ψ ∈ phiBiratAt P G A' := by
    have h2 : spanCls φ hsφ.2 ψ - (2 : ℕ) • spanCls φ hsφ.2 ψ = -spanCls φ hsφ.2 ψ := by
      rw [two_nsmul]; abel
    rw [h2] at hsub
    exact (AddSubgroup.neg_mem_iff _).mp hsub
  obtain ⟨θ, hθ⟩ := (span_iso_iff G hiso φ ψ hsφ hsψ).mpr hc
  exact ⟨θ, by rw [hθ, ← hα]⟩

/-- ★★★★**[FrdI] Theorem 5.1, (iii)** —— Frobenius-trivial な対象は **Aut-ample**。

原文 (FrdI p.99):
> and that all Frobenius-trivial objects of C are Aut-ample. In light of these
-/
theorem isAutAmple_of_frobTrivial (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {A : C} (hA : IsFrobeniusTrivial P A) : IsAutAmple P A := by
  intro g hg
  haveI := hg
  obtain ⟨θ, hθ⟩ := frobTrivial_iso_of_baseIso G hiso hA hA g
  exact ⟨θ.hom, ⟨θ.inv, θ.hom_inv_id, θ.inv_hom_id⟩, hθ⟩

/-- ★★★★**[FrdI] Theorem 5.1, (iii)** —— base-isomorphic な
Frobenius-trivial 対象は同型。

原文 (FrdI p.99):
> base-isomorphic Frobenius-trivial objects of C are, in
-/
theorem frobTrivial_baseIsomorphic_iso (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {A A' : C} (hA : IsFrobeniusTrivial P A) (hA' : IsFrobeniusTrivial P A')
    (h : BaseIsomorphic P A' A) : Nonempty (A' ≅ A) := by
  obtain ⟨e⟩ := h
  haveI : IsIso e.hom := e.isIso_hom
  obtain ⟨θ, _⟩ := frobTrivial_iso_of_baseIso G hiso hA hA' e.hom
  exact ⟨θ⟩

/-! ### ★出典の紐付け(`.src`) -/

/-- ★locator —— `Theorem 5.1, (iii)` の中心。 -/
def frobTrivial_iso_of_baseIso.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 99,
    item := "Theorem 5.1, (iii) — base-isomorphic な Frobenius-trivial 対象は同型",
    sectionId := "frdi-thm-5-1" }

end ABC3.Found.FrdI
