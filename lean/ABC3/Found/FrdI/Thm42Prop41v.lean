/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm42Prop41iv

/-!
# [FrdI] Theorem 4.2, (i) —— `Proposition 4.1, (v)` の側(コスライス)

原文 (FrdI p.77):
> paring the "primary factorizations" of x , xδ + x [cf. Deﬁnition 2.4, (i), (c), (d);

★`prop_4_1_iv` は `F` へ**入る** pre-step の側(スライス)だったが、
`prop_4_1_v` は `Do` から**出る** pre-step の側(コスライス)である。
★★コスライスでは `coaPreUnderFunctor` の値が **`P.Div` そのもの**なので、
底に沿った輸送が要らず、スライス側より簡単になる。
-/

namespace ABC3.Found.FrdI

open CategoryTheory
open scoped NNReal

universe v u w u2 v2

section Under

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-- ★★★**充満性(コスライス)** —— `Div φ ≼ Div ζ` なら `ζ` は `φ` を**通って**分解する。 -/
theorem exists_factor_under (G : Frobenioid P) {A B Z : C} (φ : A ⟶ B) (ζ : A ⟶ Z)
    (hcφ : IsCoAngular P φ) (hsφ : IsPreStep P φ)
    (hcζ : IsCoAngular P ζ) (hsζ : IsPreStep P ζ)
    (hle : MLe (P.Div φ) (P.Div ζ)) :
    ∃ ζ' : B ⟶ Z, IsCoAngular P ζ' ∧ IsPreStep P ζ' ∧ ζ = φ ≫ ζ' := by
  haveI := coaPreProp_isMultiplicative P G.core.coAngularComp
  haveI := G.coaPreUnderEquiv A
  let Zφ : Under (⟨A⟩ : WideSubcategory (coaPreProp P)) :=
    Under.mk (show (⟨A⟩ : WideSubcategory (coaPreProp P)) ⟶ ⟨B⟩ from ⟨φ, ⟨hcφ, hsφ⟩⟩)
  let Zζ : Under (⟨A⟩ : WideSubcategory (coaPreProp P)) :=
    Under.mk (show (⟨A⟩ : WideSubcategory (coaPreProp P)) ⟶ ⟨Z⟩ from ⟨ζ, ⟨hcζ, hsζ⟩⟩)
  obtain ⟨f, -⟩ := (coaPreUnderFunctor P A).map_surjective
    (show (coaPreUnderFunctor P A).obj Zφ ⟶ (coaPreUnderFunctor P A).obj Zζ from
      homOfLE (show (toOrderCat (P.Div φ) : OrderCat (Φ.val (P.toElem.obj A).base))
        ≤ toOrderCat (P.Div ζ) from hle))
  exact ⟨f.right.hom, f.right.property.1, f.right.property.2,
    (congrArg InducedWideCategory.Hom.hom (Under.w f)).symm⟩

/-- ★★★**本質的全射 ＋ 充満(コスライス)** —— `a ≼ Div φ` なら
`φ = ψ ≫ χ` と分解でき `Div ψ = a`。 -/
theorem exists_factor_of_mle_under (G : Frobenioid P) {A B : C} (φ : A ⟶ B)
    (hcφ : IsCoAngular P φ) (hsφ : IsPreStep P φ)
    {a : Φ.val (P.toElem.obj A).base} (hle : MLe a (P.Div φ)) :
    ∃ (B' : C) (ψ : A ⟶ B') (χ : B' ⟶ B), IsCoAngular P ψ ∧ IsPreStep P ψ ∧
      IsCoAngular P χ ∧ IsPreStep P χ ∧ φ = ψ ≫ χ ∧ P.Div ψ = a := by
  haveI := coaPreProp_isMultiplicative P G.core.coAngularComp
  haveI := G.coaPreUnderEquiv A
  let Zφ : Under (⟨A⟩ : WideSubcategory (coaPreProp P)) :=
    Under.mk (show (⟨A⟩ : WideSubcategory (coaPreProp P)) ⟶ ⟨B⟩ from ⟨φ, ⟨hcφ, hsφ⟩⟩)
  let W := (coaPreUnderFunctor P A).objPreimage (toOrderCat a)
  have hiso : (coaPreUnderFunctor P A).obj W ≅ toOrderCat a :=
    (coaPreUnderFunctor P A).objObjPreimageIso _
  have hWa : P.Div W.hom.hom = a :=
    mle_antisymm (P.divisorial _).1.1 (P.divisorial _).2
      (leOfHom hiso.hom) (leOfHom hiso.inv)
  obtain ⟨f, -⟩ := (coaPreUnderFunctor P A).map_surjective
    (show (coaPreUnderFunctor P A).obj W ⟶ (coaPreUnderFunctor P A).obj Zφ from
      hiso.hom ≫ homOfLE (show (toOrderCat a : OrderCat (Φ.val (P.toElem.obj A).base))
        ≤ toOrderCat (P.Div φ) from hle))
  exact ⟨W.right.obj, W.hom.hom, f.right.hom, W.hom.property.1, W.hom.property.2,
    f.right.property.1, f.right.property.2,
    (congrArg InducedWideCategory.Hom.hom (Under.w f)).symm, hWa⟩

/-! ## ★★★★★コスライス版の「同素性」 -/

/-- ★★**「同じ素点」の圏論的条件(コスライス版)**。 -/
def SamePrimeCatOut (P : PreFrobenioid C Φ) {A B B' : C} (ϵ : A ⟶ B) (ϵ' : A ⟶ B') : Prop :=
  ∃ (Z : C) (ζ : A ⟶ Z), IsCoAngular P ζ ∧ IsStep P ζ ∧
    (∃ ψ : Z ⟶ B, IsPreStep P ψ ∧ ϵ = ζ ≫ ψ) ∧
    (∃ ψ' : Z ⟶ B', IsPreStep P ψ' ∧ ϵ' = ζ ≫ ψ')

/-- ★`IsStep` は `Div ≠ 0` と同じ(isotropic 型)。 -/
theorem isStep_iff_div_ne_zero (hiso : ∀ X : C, IsIsotropic P X)
    {A B : C} (φ : B ⟶ A) (h : IsPreStep P φ) : IsStep P φ ↔ P.Div φ ≠ 0 := by
  rw [isStep_iff_preStepVal_ne_zero hiso φ h, Ne, Ne, preStepVal_eq_zero_iff]

/-- ★★★★★★**コスライス版の同素性判定**。 -/
theorem toPrime_div_eq_iff (G : Frobenioid P) (hiso : ∀ X : C, IsIsotropic P X) {A : C}
    {ι₀ : Prime (Φ.val (P.toElem.obj A).base) → Pf (Φ.val (P.toElem.obj A).base) → ℝ≥0}
    (H : IsPerfFactorialWith (Φ.val (P.toElem.obj A).base) ι₀)
    (hperf : IsPerfectMonoid (Φ.val (P.toElem.obj A).base))
    (hdivM : IsDivisorial (Φ.val (P.toElem.obj A).base))
    {B B' : C} (ϵ : A ⟶ B) (hcϵ : IsCoAngular P ϵ) (hsϵ : IsPreStep P ϵ)
    (hpϵ : IsPrimaryElt (P.Div ϵ))
    (ϵ' : A ⟶ B') (hcϵ' : IsCoAngular P ϵ') (hsϵ' : IsPreStep P ϵ')
    (hpϵ' : IsPrimaryElt (P.Div ϵ')) :
    toPrime _ (P.Div ϵ) hpϵ = toPrime _ (P.Div ϵ') hpϵ' ↔ SamePrimeCatOut P ϵ ϵ' := by
  have hsup : SuppElt ι₀ (P.Div ϵ) = SuppElt ι₀ (P.Div ϵ')
      ↔ toPrime _ (P.Div ϵ) hpϵ = toPrime _ (P.Div ϵ') hpϵ' := by
    rw [suppElt_eq_singleton_toPrime H hperf hdivM hpϵ,
      suppElt_eq_singleton_toPrime H hperf hdivM hpϵ']
    exact ⟨fun h => Set.singleton_injective h, fun h => by rw [h]⟩
  rw [← hsup, suppElt_eq_iff_exists_common H hperf hdivM hpϵ hpϵ']
  constructor
  · rintro ⟨d, hd0, hdϵ, hdϵ'⟩
    obtain ⟨Z, ζ, χ, hζc, hζs, hχc, hχs, hfac, hζval⟩ :=
      exists_factor_of_mle_under G ϵ hcϵ hsϵ hdϵ
    have hle' : MLe (P.Div ζ) (P.Div ϵ') := by rw [hζval]; exact hdϵ'
    obtain ⟨ψ', hψ'c, hψ's, hfac'⟩ :=
      exists_factor_under G ζ ϵ' hζc hζs hcϵ' hsϵ' hle'
    refine ⟨Z, ζ, hζc, ?_, ⟨χ, hχs, hfac⟩, ⟨ψ', hψ's, hfac'⟩⟩
    refine (isStep_iff_div_ne_zero hiso ζ hζs).mpr ?_
    rw [hζval]; exact hd0
  · rintro ⟨Z, ζ, hζc, hζst, ⟨ψ, hψs, hfac⟩, ⟨ψ', hψ's, hfac'⟩⟩
    refine ⟨P.Div ζ, (isStep_iff_div_ne_zero hiso ζ hζst.1).mp hζst, ?_, ?_⟩
    · rw [hfac]; exact mle_div_comp ζ ψ hψs.1
    · rw [hfac']; exact mle_div_comp ζ ψ' hψ's.1

def SamePrimeCatOut.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 80,
    item := "Theorem 4.2, (ii) — 同素性の圏論的判定（コスライス）",
    sectionId := "frdi-thm-4-2" }

/-! ## ★★★★★★`Proposition 4.1, (v)` の圏論版 -/

/-- ★★どの素点も `𝒪^▷` の元の零因子で実現される。 -/
theorem exists_otri_realizing_div
    (hdivS : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
      ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
    {A : C} (p : Prime (Φ.val (P.toElem.obj A).base)) :
    ∃ (u : End A) (hu : u ∈ OTri P A) (hp : IsPrimaryElt (P.Div (((u : A ⟶ A))))),
      toPrime _ (P.Div (((u : A ⟶ A)))) hp = p := by
  obtain ⟨⟨a, ha⟩, hx⟩ := Quotient.exists_rep p
  obtain ⟨u, hu⟩ := hdivS A a
  subst hu
  exact ⟨(u : End A), u.2, ha, hx⟩

/-- ★★★★**`Proposition 4.1, (v)` の右辺の圏論版**。 -/
def Prop41vRhsCat (P : PreFrobenioid C Φ) {Do E F : C} (δ : Do ⟶ E) (ϵ : E ⟶ F) : Prop :=
  ∃ (E₀ : C) (δ₀ : Do ⟶ E₀), IsCoAngular P δ₀ ∧ IsPreStep P δ₀ ∧ IsPrimaryElt (P.Div δ₀) ∧
    ∀ (E' : C) (δ' : Do ⟶ E'), IsCoAngular P δ' → IsPreStep P δ' → IsPrimaryElt (P.Div δ') →
      ¬ SamePrimeCatOut P δ' δ₀ →
        ((∃ ζ : E' ⟶ E, IsPreStep P ζ ∧ δ = δ' ≫ ζ) ↔
         (∃ θ : E' ⟶ F, IsPreStep P θ ∧ δ ≫ ϵ = δ' ≫ θ))

/-- ★★★★★★**[FrdI] Proposition 4.1, (v) の圏論版**。 -/
theorem prop_4_1_v_cat (G : Frobenioid P) (hiso : ∀ X : C, IsIsotropic P X)
    (hdivS : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
      ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
    {Do E F : C} (δ : Do ⟶ E) (ϵ : E ⟶ F)
    (hsδ : IsPreStep P δ) (hsϵ : IsPreStep P ϵ) (hϵ0 : P.Div ϵ ≠ 0)
    {ι₀ : Prime (Φ.val (P.toElem.obj Do).base) → Pf (Φ.val (P.toElem.obj Do).base) → ℝ≥0}
    (H : IsPerfFactorialWith (Φ.val (P.toElem.obj Do).base) ι₀)
    (hperf : IsPerfectMonoid (Φ.val (P.toElem.obj Do).base))
    (hdivM : IsDivisorial (Φ.val (P.toElem.obj Do).base)) :
    IsPrimaryElt (P.Div ϵ) ↔ Prop41vRhsCat P δ ϵ := by
  rw [prop_4_1_v G hiso δ ϵ hsδ hsϵ hϵ0 H hperf]
  constructor
  · rintro ⟨p, hp⟩
    obtain ⟨u, hu, hup, hutp⟩ := exists_otri_realizing_div hdivS p
    refine ⟨Do, ((u : End Do) : Do ⟶ Do), prop_1_4_i P _ (fun Y _ => hiso Y),
      isPreStep_of_otri _ hu, hup, ?_⟩
    intro E' δ' hcδ' hsδ' hpδ' hns
    refine hp E' δ' hsδ' hpδ' ?_
    intro hsupp
    refine hns ?_
    refine (toPrime_div_eq_iff G hiso H hperf hdivM δ' hcδ' hsδ' hpδ'
      ((u : End Do) : Do ⟶ Do) (prop_1_4_i P _ (fun Y _ => hiso Y))
      (isPreStep_of_otri _ hu) hup).mp ?_
    rw [suppElt_eq_singleton_toPrime H hperf hdivM hpδ'] at hsupp
    have h2 : toPrime _ (P.Div δ') hpδ' = p := Set.singleton_injective hsupp
    rw [h2, ← hutp]
  · rintro ⟨E₀, δ₀, hc₀, hs₀, hp₀, hcond⟩
    refine ⟨toPrime _ (P.Div δ₀) hp₀, ?_⟩
    intro E' δ' hsδ' hpδ' hne
    have hcδ' : IsCoAngular P δ' := prop_1_4_i P _ (fun Y _ => hiso Y)
    refine hcond E' δ' hcδ' hsδ' hpδ' ?_
    intro hsame
    refine hne ?_
    have h3 : toPrime _ (P.Div δ') hpδ' = toPrime _ (P.Div δ₀) hp₀ :=
      (toPrime_div_eq_iff G hiso H hperf hdivM δ' hcδ' hsδ' hpδ' δ₀ hc₀ hs₀ hp₀).mpr hsame
    rw [suppElt_eq_singleton_toPrime H hperf hdivM hpδ', h3]

def prop_4_1_v_cat.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 77,
    item := "Proposition 4.1, (v) — 右辺の圏論版",
    sectionId := "frdi-prop-4-1" }

end Under

end ABC3.Found.FrdI
