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

/-! ## ★★★★★★「入る primary」と「出る primary」の橋

★★`Div : 𝒪^▷ ↠ Φ` により、`A` の零因子はすべて **`A` の自己射**で実現される。
★自己射は「`A` へ入る pre-step」でも「`A` から出る pre-step」でもあるので、
これが 2 つの側をつなぐ。 -/

section Bridge

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}

/-- ★同型を後置しても `Div` は変わらない(pre-step のとき)。 -/
theorem div_comp_iso {X Y Z : C₂} (f : X ⟶ Y) (e : Y ⟶ Z) [IsIso e] :
    P₂.Div (f ≫ e) = P₂.Div f := by
  rw [P₂.Div_comp, show P₂.Div e = 0 from isIsometric_of_isIso P₂ e, map_zero, zero_add,
    show P₂.degFr e = 1 from degFr_of_isIso P₂ e]
  exact one_smul _ _

/-- ★★★★★**「`A` へ入る primary の保存」⟹「`A` から出る primary の保存」**。 -/
theorem primaryOut_of_primaryIn (Ψ : C ⥤ C₂) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hdivS : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
      ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
    {A : C}
    (hIn : ∀ (E' : C) (ϵ' : E' ⟶ A), IsPreStep P ϵ' → IsPrimaryElt (P.Div ϵ') →
      IsPrimaryElt (P₂.Div (Ψ.map ϵ')))
    {B : C} (δ : A ⟶ B) (hsδ : IsPreStep P δ) (hpδ : IsPrimaryElt (P.Div δ)) :
    IsPrimaryElt (P₂.Div (Ψ.map δ)) := by
  obtain ⟨u, hu⟩ := hdivS A (P.Div δ)
  have hcu : IsCoAngular P (((u : End A)) : A ⟶ A) := prop_1_4_i P _ (fun Y _ => hiso Y)
  have hsu : IsPreStep P (((u : End A)) : A ⟶ A) := isPreStep_of_otri _ u.2
  have hcδ : IsCoAngular P δ := prop_1_4_i P _ (fun Y _ => hiso Y)
  obtain ⟨θ, hθiso, hθ⟩ :=
    exists_iso_of_div_eq G (((u : End A)) : A ⟶ A) δ hcu hsu hcδ hsδ hu
  haveI := hθiso
  have hpu : IsPrimaryElt (P.Div (((u : End A)) : A ⟶ A)) := by rw [hu]; exact hpδ
  have h1 : IsPrimaryElt (P₂.Div (Ψ.map (((u : End A)) : A ⟶ A))) := hIn A _ hsu hpu
  have h2 : P₂.Div (Ψ.map δ) = P₂.Div (Ψ.map (((u : End A)) : A ⟶ A)) := by
    rw [← hθ, Ψ.map_comp]
    haveI : IsIso (Ψ.map θ) := Ψ.map_isIso θ
    exact div_comp_iso (Ψ.map (((u : End A)) : A ⟶ A)) (Ψ.map θ)
  rw [h2]
  exact h1

def primaryOut_of_primaryIn.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 79,
    item := "Theorem 4.2, (i) — 入る primary と出る primary の橋",
    sectionId := "frdi-thm-4-2" }

/-- ★★★**「後置の値」が等しい 2 本の co-angular pre-step は始域の上で同型**。

★`exists_iso_of_div_eq` のスライス版(`Definition 1.3, (iii), (d)` の**後置**の圏同値)。 -/
theorem exists_iso_of_preStepVal_eq (G : Frobenioid P) {A B B' : C} (φ : B ⟶ A) (φ' : B' ⟶ A)
    (hcφ : IsCoAngular P φ) (hsφ : IsPreStep P φ)
    (hcφ' : IsCoAngular P φ') (hsφ' : IsPreStep P φ')
    (h : preStepVal P φ hsφ = preStepVal P φ' hsφ') :
    ∃ θ : B ⟶ B', IsIso θ ∧ θ ≫ φ' = φ := by
  haveI := coaPreProp_isMultiplicative P G.core.coAngularComp
  haveI := G.coaPreOverEquiv A
  let Z : Over (⟨A⟩ : WideSubcategory (coaPreProp P)) :=
    Over.mk (show (⟨B⟩ : WideSubcategory (coaPreProp P)) ⟶ ⟨A⟩ from ⟨φ, hcφ, hsφ⟩)
  let W : Over (⟨A⟩ : WideSubcategory (coaPreProp P)) :=
    Over.mk (show (⟨B'⟩ : WideSubcategory (coaPreProp P)) ⟶ ⟨A⟩ from ⟨φ', hcφ', hsφ'⟩)
  have hobj : (coaPreOverFunctor P A).obj Z = (coaPreOverFunctor P A).obj W := by
    show Opposite.op (toOrderCat (preStepVal P φ hsφ))
      = Opposite.op (toOrderCat (preStepVal P φ' hsφ'))
    rw [h]
  obtain ⟨f, -⟩ := (coaPreOverFunctor P A).map_surjective (eqToHom hobj)
  obtain ⟨g, -⟩ := (coaPreOverFunctor P A).map_surjective (eqToHom hobj.symm)
  have hfg : f ≫ g = 𝟙 Z :=
    (coaPreOverFunctor P A).map_injective (Subsingleton.elim _ _)
  have hgf : g ≫ f = 𝟙 W :=
    (coaPreOverFunctor P A).map_injective (Subsingleton.elim _ _)
  refine ⟨f.left.hom, ⟨g.left.hom, ?_, ?_⟩, ?_⟩
  · exact congrArg (fun t : Z ⟶ Z => t.left.hom) hfg
  · exact congrArg (fun t : W ⟶ W => t.left.hom) hgf
  · exact congrArg InducedWideCategory.Hom.hom (Over.w f)

/-- ★★★★★**「`A` から出る primary の保存」⟹「`A` へ入る primary の保存」**。 -/
theorem primaryIn_of_primaryOut (Ψ : C ⥤ C₂) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hdivS : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
      ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
    (hPS : ∀ {X Y : C} (f : X ⟶ Y), IsPreStep P f → IsPreStep P₂ (Ψ.map f))
    {A : C}
    (hOut : ∀ (B' : C) (δ' : A ⟶ B'), IsPreStep P δ' → IsPrimaryElt (P.Div δ') →
      IsPrimaryElt (P₂.Div (Ψ.map δ')))
    {E : C} (ϵ : E ⟶ A) (hsϵ : IsPreStep P ϵ) (hpϵ : IsPrimaryElt (P.Div ϵ)) :
    IsPrimaryElt (P₂.Div (Ψ.map ϵ)) := by
  obtain ⟨u, hu⟩ := hdivS A (preStepVal P ϵ hsϵ)
  have hsu : IsPreStep P (((u : End A)) : A ⟶ A) := isPreStep_of_otri _ u.2
  have hcu : IsCoAngular P (((u : End A)) : A ⟶ A) := prop_1_4_i P _ (fun Y _ => hiso Y)
  have hcϵ : IsCoAngular P ϵ := prop_1_4_i P _ (fun Y _ => hiso Y)
  have hvu : preStepVal P (((u : End A)) : A ⟶ A) hsu = preStepVal P ϵ hsϵ := by
    rw [preStepVal_of_otri _ u.2]; exact hu
  obtain ⟨θ, hθiso, hθ⟩ :=
    exists_iso_of_preStepVal_eq G (((u : End A)) : A ⟶ A) ϵ hcu hsu hcϵ hsϵ hvu
  haveI := hθiso
  have hpu : IsPrimaryElt (P.Div (((u : End A)) : A ⟶ A)) := by
    rw [hu]; exact (isPrimaryElt_preStepVal_iff ϵ hsϵ).mpr hpϵ
  have h1 : IsPrimaryElt (P₂.Div (Ψ.map (((u : End A)) : A ⟶ A))) := hOut A _ hsu hpu
  rw [← hθ, Ψ.map_comp] at h1
  haveI : IsIso (Ψ.map θ) := Ψ.map_isIso θ
  exact (isPrimaryElt_div_iso_comp (Ψ.map θ) (Ψ.map ϵ) (hPS ϵ hsϵ)).mp h1

def primaryIn_of_primaryOut.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 79,
    item := "Theorem 4.2, (i) — 出る primary と入る primary の橋",
    sectionId := "frdi-thm-4-2" }

/-! ## ★★★★★v 側の移送 -/

/-- ★`f` が `g` から出て pre-step で分解する。 -/
def FactorsPost (P : PreFrobenioid C Φ) {X Y Z : C} (f : X ⟶ Z) (g : X ⟶ Y) : Prop :=
  ∃ ζ : Y ⟶ Z, IsPreStep P ζ ∧ f = g ≫ ζ

/-- ★同型を後置しても変わらない。 -/
theorem factorsPost_iso_comp {X Y Y' Z : C} (e : Y' ⟶ Y) [IsIso e] (f : X ⟶ Z) (g : X ⟶ Y') :
    FactorsPost P f (g ≫ e) ↔ FactorsPost P f g := by
  constructor
  · rintro ⟨ζ, hζ, hfac⟩
    exact ⟨e ≫ ζ, IsPreStep.comp P (isPreStep_of_isIso P e) hζ, by rw [hfac, Category.assoc]⟩
  · rintro ⟨ζ, hζ, hfac⟩
    refine ⟨inv e ≫ ζ, IsPreStep.comp P (isPreStep_of_isIso P (inv e)) hζ, ?_⟩
    rw [hfac, Category.assoc, ← Category.assoc e, IsIso.hom_inv_id, Category.id_comp]

/-- ★★「pre-step から出る分解」は圏同値で行き来する。 -/
theorem factorsPost_map_iff (Ψ : C ≌ C₂)
    (hPS : ∀ {X Y : C} (f : X ⟶ Y), IsPreStep P f → IsPreStep P₂ (Ψ.functor.map f))
    (hPS' : ∀ {X Y : C} (f : X ⟶ Y), IsPreStep P₂ (Ψ.functor.map f) → IsPreStep P f)
    {X Y Z : C} (f : X ⟶ Z) (g : X ⟶ Y) :
    FactorsPost P f g ↔ FactorsPost P₂ (Ψ.functor.map f) (Ψ.functor.map g) := by
  constructor
  · rintro ⟨ζ, hζ, hfac⟩
    exact ⟨Ψ.functor.map ζ, hPS ζ hζ, by rw [hfac, Ψ.functor.map_comp]⟩
  · rintro ⟨ζ', hζ', hfac⟩
    obtain ⟨ζ, hζeq⟩ := Ψ.functor.map_surjective ζ'
    refine ⟨ζ, hPS' ζ (by rw [hζeq]; exact hζ'), ?_⟩
    refine Ψ.functor.map_injective ?_
    rw [Ψ.functor.map_comp, hζeq, hfac]

/-- ★★★★★**`SamePrimeCatOut` は `Ψ` で移る**。 -/
theorem samePrimeCatOut_map (Ψ : C ⥤ C₂)
    (hCA : ∀ {X Y : C} (f : X ⟶ Y), IsCoAngular P f → IsCoAngular P₂ (Ψ.map f))
    (hPS : ∀ {X Y : C} (f : X ⟶ Y), IsPreStep P f → IsPreStep P₂ (Ψ.map f))
    (hST : ∀ {X Y : C} (f : X ⟶ Y), IsStep P f → IsStep P₂ (Ψ.map f))
    {A B B' : C} {ϵ : A ⟶ B} {ϵ' : A ⟶ B'} (h : SamePrimeCatOut P ϵ ϵ') :
    SamePrimeCatOut P₂ (Ψ.map ϵ) (Ψ.map ϵ') := by
  obtain ⟨Z, ζ, hζc, hζst, ⟨ψ, hψs, hfac⟩, ⟨ψ', hψ's, hfac'⟩⟩ := h
  exact ⟨Ψ.obj Z, Ψ.map ζ, hCA ζ hζc, hST ζ hζst,
    ⟨Ψ.map ψ, hPS ψ hψs, by rw [hfac, Ψ.map_comp]⟩,
    ⟨Ψ.map ψ', hPS ψ' hψ's, by rw [hfac', Ψ.map_comp]⟩⟩

/-- ★★同型を後置しても `SamePrimeCatOut` は変わらない。 -/
theorem samePrimeCatOut_iso_comp {A B B' X : C} (e : B ⟶ X) [IsIso e]
    {ϵ : A ⟶ B} {ϵ' : A ⟶ B'} :
    SamePrimeCatOut P (ϵ ≫ e) ϵ' ↔ SamePrimeCatOut P ϵ ϵ' := by
  constructor
  · rintro ⟨Z, ζ, hζc, hζst, ⟨ψ, hψs, hfac⟩, hr⟩
    refine ⟨Z, ζ, hζc, hζst, ⟨ψ ≫ inv e, IsPreStep.comp P hψs (isPreStep_of_isIso P (inv e)),
      ?_⟩, hr⟩
    rw [← Category.assoc, ← hfac, Category.assoc, IsIso.hom_inv_id, Category.comp_id]
  · rintro ⟨Z, ζ, hζc, hζst, ⟨ψ, hψs, hfac⟩, hr⟩
    refine ⟨Z, ζ, hζc, hζst, ⟨ψ ≫ e, IsPreStep.comp P hψs (isPreStep_of_isIso P e), ?_⟩, hr⟩
    rw [← Category.assoc, ← hfac]

/-- ★★★★★★**`Proposition 4.1, (v)` の右辺は `Ψ` で移る**。 -/
theorem prop41vRhsCat_map (Ψ : C ≌ C₂)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hCA : ∀ {X Y : C} (f : X ⟶ Y), IsCoAngular P f → IsCoAngular P₂ (Ψ.functor.map f))
    (hPS : ∀ {X Y : C} (f : X ⟶ Y), IsPreStep P f → IsPreStep P₂ (Ψ.functor.map f))
    (hPS' : ∀ {X Y : C} (f : X ⟶ Y), IsPreStep P₂ (Ψ.functor.map f) → IsPreStep P f)
    (hST : ∀ {X Y : C} (f : X ⟶ Y), IsStep P f → IsStep P₂ (Ψ.functor.map f))
    {Do E F : C} {δ : Do ⟶ E} {ϵ : E ⟶ F}
    (hFwd : ∀ (E' : C) (δ' : Do ⟶ E'), IsPreStep P δ' →
      IsPrimaryElt (P.Div δ') → IsPrimaryElt (P₂.Div (Ψ.functor.map δ')))
    (hBwd : ∀ (E' : C) (δ' : Do ⟶ E'), IsPreStep P δ' →
      IsPrimaryElt (P₂.Div (Ψ.functor.map δ')) → IsPrimaryElt (P.Div δ'))
    (h : Prop41vRhsCat P δ ϵ) :
    Prop41vRhsCat P₂ (Ψ.functor.map δ) (Ψ.functor.map ϵ) := by
  obtain ⟨E₀, δ₀, hc₀, hs₀, hp₀, hcond⟩ := h
  refine ⟨Ψ.functor.obj E₀, Ψ.functor.map δ₀, hCA δ₀ hc₀, hPS δ₀ hs₀,
    hFwd E₀ δ₀ hs₀ hp₀, ?_⟩
  intro E' δ' hcδ' hsδ' hpδ' hns
  obtain ⟨E'₀, ⟨e⟩⟩ := Functor.EssSurj.mem_essImage (F := Ψ.functor) E'
  obtain ⟨δ'₀, hδ'₀⟩ := Ψ.functor.map_surjective (δ' ≫ e.inv)
  have hsδ'₀ : IsPreStep P δ'₀ := by
    refine hPS' δ'₀ ?_
    rw [hδ'₀]
    exact IsPreStep.comp P₂ hsδ' (isPreStep_of_isIso P₂ e.inv)
  have hpδ'₀ : IsPrimaryElt (P.Div δ'₀) := by
    refine hBwd E'₀ δ'₀ hsδ'₀ ?_
    rw [hδ'₀]
    haveI : IsIso e.inv := inferInstance
    rw [div_comp_iso δ' e.inv]
    exact hpδ'
  have hns₀ : ¬ SamePrimeCatOut P δ'₀ δ₀ := by
    intro hc
    refine hns ?_
    have h1 := samePrimeCatOut_map Ψ.functor hCA hPS hST hc
    rw [hδ'₀] at h1
    exact (samePrimeCatOut_iso_comp e.inv).mp h1
  have hiff := hcond E'₀ δ'₀ (prop_1_4_i P _ (fun Y _ => hiso Y)) hsδ'₀ hpδ'₀ hns₀
  have hL : FactorsPost P₂ (Ψ.functor.map δ) δ' ↔ FactorsPost P δ δ'₀ := by
    rw [factorsPost_map_iff Ψ hPS hPS' δ δ'₀, hδ'₀]
    exact (factorsPost_iso_comp e.inv (Ψ.functor.map δ) δ').symm
  have hR : FactorsPost P₂ (Ψ.functor.map δ ≫ Ψ.functor.map ϵ) δ'
      ↔ FactorsPost P (δ ≫ ϵ) δ'₀ := by
    rw [factorsPost_map_iff Ψ hPS hPS' (δ ≫ ϵ) δ'₀, hδ'₀, Ψ.functor.map_comp]
    exact (factorsPost_iso_comp e.inv (Ψ.functor.map δ ≫ Ψ.functor.map ϵ) δ').symm
  exact hL.trans (hiff.trans hR.symm)

def prop41vRhsCat_map.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 79,
    item := "Theorem 4.2, (i) — Proposition 4.1, (v) の右辺は Ψ で移る",
    sectionId := "frdi-thm-4-2" }

/-- ★★★★★★**primary 性の「一段の拡張」(v 側)** ——
`Do` から出る pre-step について primary 性が両向きに移るなら、
`Do` から pre-step `δ : Do ⟶ E` を持つ `E` から出る pre-step についても移る。 -/
theorem isPrimaryElt_div_extend_v (Ψ : C ≌ C₂)
    (G : Frobenioid P) (G₂ : Frobenioid P₂)
    (hiso : ∀ X : C, IsIsotropic P X) (hiso₂ : ∀ X : C₂, IsIsotropic P₂ X)
    (hdivS : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
      ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
    (hdivS₂ : ∀ (Y : C₂) (a : Φ₂.val (P₂.toElem.obj Y).base),
      ∃ u : OTri P₂ Y, P₂.Div (((u : End Y) : Y ⟶ Y)) = a)
    (hCA : ∀ {X Y : C} (f : X ⟶ Y), IsCoAngular P f → IsCoAngular P₂ (Ψ.functor.map f))
    (hPS : ∀ {X Y : C} (f : X ⟶ Y), IsPreStep P f → IsPreStep P₂ (Ψ.functor.map f))
    (hPS' : ∀ {X Y : C} (f : X ⟶ Y), IsPreStep P₂ (Ψ.functor.map f) → IsPreStep P f)
    (hST : ∀ {X Y : C} (f : X ⟶ Y), IsStep P f → IsStep P₂ (Ψ.functor.map f))
    {Do E F : C} (δ : Do ⟶ E) (ϵ : E ⟶ F)
    (hsδ : IsPreStep P δ) (hsϵ : IsPreStep P ϵ)
    (hϵ0 : P.Div ϵ ≠ 0) (hϵ0₂ : P₂.Div (Ψ.functor.map ϵ) ≠ 0)
    {ι₀ : Prime (Φ.val (P.toElem.obj Do).base) → Pf (Φ.val (P.toElem.obj Do).base) → ℝ≥0}
    (H : IsPerfFactorialWith (Φ.val (P.toElem.obj Do).base) ι₀)
    (hperf : IsPerfectMonoid (Φ.val (P.toElem.obj Do).base))
    (hdivM : IsDivisorial (Φ.val (P.toElem.obj Do).base))
    {ι₂ : Prime (Φ₂.val (P₂.toElem.obj (Ψ.functor.obj Do)).base) →
      Pf (Φ₂.val (P₂.toElem.obj (Ψ.functor.obj Do)).base) → ℝ≥0}
    (H₂ : IsPerfFactorialWith (Φ₂.val (P₂.toElem.obj (Ψ.functor.obj Do)).base) ι₂)
    (hperf₂ : IsPerfectMonoid (Φ₂.val (P₂.toElem.obj (Ψ.functor.obj Do)).base))
    (hdivM₂ : IsDivisorial (Φ₂.val (P₂.toElem.obj (Ψ.functor.obj Do)).base))
    (hFwd : ∀ (E' : C) (δ' : Do ⟶ E'), IsPreStep P δ' →
      IsPrimaryElt (P.Div δ') → IsPrimaryElt (P₂.Div (Ψ.functor.map δ')))
    (hBwd : ∀ (E' : C) (δ' : Do ⟶ E'), IsPreStep P δ' →
      IsPrimaryElt (P₂.Div (Ψ.functor.map δ')) → IsPrimaryElt (P.Div δ'))
    (h : IsPrimaryElt (P.Div ϵ)) :
    IsPrimaryElt (P₂.Div (Ψ.functor.map ϵ)) :=
  (prop_4_1_v_cat G₂ hiso₂ hdivS₂ (Ψ.functor.map δ) (Ψ.functor.map ϵ)
      (hPS δ hsδ) (hPS ϵ hsϵ) hϵ0₂ H₂ hperf₂ hdivM₂).mpr
    (prop41vRhsCat_map Ψ hiso hCA hPS hPS' hST hFwd hBwd
      ((prop_4_1_v_cat G hiso hdivS δ ϵ hsδ hsϵ hϵ0 H hperf hdivM).mp h))

def isPrimaryElt_div_extend_v.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 79,
    item := "Theorem 4.2, (i) — primary 性の一段の拡張（v 側）",
    sectionId := "frdi-thm-4-2" }

/-! ## ★★★★★★橋の **↔** 版(保存かつ反射) -/

/-- ★★★★★★**「`A` へ入る primary が両向きに移る」⟹「`A` から出る primary も両向きに移る」**。 -/
theorem primaryOut_iff_of_primaryIn (Ψ : C ⥤ C₂) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hdivS : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
      ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
    {A : C}
    (hIn : ∀ (E' : C) (ϵ' : E' ⟶ A), IsPreStep P ϵ' →
      (IsPrimaryElt (P.Div ϵ') ↔ IsPrimaryElt (P₂.Div (Ψ.map ϵ'))))
    {B : C} (δ : A ⟶ B) (hsδ : IsPreStep P δ) :
    IsPrimaryElt (P.Div δ) ↔ IsPrimaryElt (P₂.Div (Ψ.map δ)) := by
  obtain ⟨u, hu⟩ := hdivS A (P.Div δ)
  have hcu : IsCoAngular P (((u : End A)) : A ⟶ A) := prop_1_4_i P _ (fun Y _ => hiso Y)
  have hsu : IsPreStep P (((u : End A)) : A ⟶ A) := isPreStep_of_otri _ u.2
  have hcδ : IsCoAngular P δ := prop_1_4_i P _ (fun Y _ => hiso Y)
  obtain ⟨θ, hθiso, hθ⟩ :=
    exists_iso_of_div_eq G (((u : End A)) : A ⟶ A) δ hcu hsu hcδ hsδ hu
  haveI := hθiso
  have h2 : P₂.Div (Ψ.map δ) = P₂.Div (Ψ.map (((u : End A)) : A ⟶ A)) := by
    rw [← hθ, Ψ.map_comp]
    haveI : IsIso (Ψ.map θ) := Ψ.map_isIso θ
    exact div_comp_iso (Ψ.map (((u : End A)) : A ⟶ A)) (Ψ.map θ)
  rw [h2, ← hu]
  exact hIn A _ hsu

/-- ★★★★★★**「`A` から出る primary が両向きに移る」⟹「`A` へ入る primary も両向きに移る」**。 -/
theorem primaryIn_iff_of_primaryOut (Ψ : C ⥤ C₂) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hdivS : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
      ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
    (hPS : ∀ {X Y : C} (f : X ⟶ Y), IsPreStep P f → IsPreStep P₂ (Ψ.map f))
    {A : C}
    (hOut : ∀ (B' : C) (δ' : A ⟶ B'), IsPreStep P δ' →
      (IsPrimaryElt (P.Div δ') ↔ IsPrimaryElt (P₂.Div (Ψ.map δ'))))
    {E : C} (ϵ : E ⟶ A) (hsϵ : IsPreStep P ϵ) :
    IsPrimaryElt (P.Div ϵ) ↔ IsPrimaryElt (P₂.Div (Ψ.map ϵ)) := by
  obtain ⟨u, hu⟩ := hdivS A (preStepVal P ϵ hsϵ)
  have hsu : IsPreStep P (((u : End A)) : A ⟶ A) := isPreStep_of_otri _ u.2
  have hcu : IsCoAngular P (((u : End A)) : A ⟶ A) := prop_1_4_i P _ (fun Y _ => hiso Y)
  have hcϵ : IsCoAngular P ϵ := prop_1_4_i P _ (fun Y _ => hiso Y)
  have hvu : preStepVal P (((u : End A)) : A ⟶ A) hsu = preStepVal P ϵ hsϵ := by
    rw [preStepVal_of_otri _ u.2]; exact hu
  obtain ⟨θ, hθiso, hθ⟩ :=
    exists_iso_of_preStepVal_eq G (((u : End A)) : A ⟶ A) ϵ hcu hsu hcϵ hsϵ hvu
  haveI := hθiso
  have hleft : IsPrimaryElt (P.Div ϵ)
      ↔ IsPrimaryElt (P.Div (((u : End A)) : A ⟶ A)) := by
    rw [hu]
    exact (isPrimaryElt_preStepVal_iff ϵ hsϵ).symm
  have hright : IsPrimaryElt (P₂.Div (Ψ.map (((u : End A)) : A ⟶ A)))
      ↔ IsPrimaryElt (P₂.Div (Ψ.map ϵ)) := by
    rw [← hθ, Ψ.map_comp]
    haveI : IsIso (Ψ.map θ) := Ψ.map_isIso θ
    exact isPrimaryElt_div_iso_comp (Ψ.map θ) (Ψ.map ϵ) (hPS ϵ hsϵ)
  exact hleft.trans ((hOut A _ hsu).trans hright)

def primaryOut_iff_of_primaryIn.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 79,
    item := "Theorem 4.2, (i) — 入る/出る primary の橋（両向き）",
    sectionId := "frdi-thm-4-2" }

end Bridge

end ABC3.Found.FrdI
