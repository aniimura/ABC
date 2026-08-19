/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm42PrimeCat

/-!
# [FrdI] Theorem 4.2, (i) —— `Proposition 4.1, (iv)` の右辺を**圏論的**にする

原文 (FrdI p.77):
> For every primary element x ∈/ p, x ≤ x if and only if x ≤ xδ + x .

★★在庫の `prop_4_1_iv` の右辺には `Prime` と `SuppElt` が現れる:
```
∃ p : Prime (Φ(F)), ∀ E' ϵ' , IsPrimaryElt (val ϵ') → SuppElt (val ϵ') ≠ {p} → (…)
```
★これが「圏論的でない」ため、`Ψ` へそのまま移らなかった。

## ★★★★★★2 つの在庫でこれが解ける

1. `suppElt_eq_singleton_toPrime`(在庫) —— `SuppElt (val ϵ') = {toPrime (val ϵ')}`
   したがって `SuppElt (val ϵ') ≠ {p}` は **`toPrime (val ϵ') ≠ p`** と同じ。
2. `toPrime_preStepVal_eq_iff`(`Thm42PrimeCat.lean`) ——
   `toPrime` の等式は **`SamePrimeCat`**(圏論的)と同値。

★★さらに `Div : 𝒪^▷ ↠ Φ` により **どの素点も primary pre-step で実現される**ので、
`∃ p : Prime (Φ(F))` は `∃ ϵ₀ : primary な pre-step into F` に書き換わる。

★★★結論: **`Proposition 4.1, (iv)` の右辺は圏論的に書ける**。
-/

namespace ABC3.Found.FrdI

open CategoryTheory
open scoped NNReal

universe v u w u2 v2

section Cat41iv

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-- ★★★★**`Proposition 4.1, (iv)` の右辺の圏論版**。 -/
def Prop41ivRhsCat (P : PreFrobenioid C Φ) {Do E F : C} (δ : Do ⟶ E) (ϵ : E ⟶ F) : Prop :=
  ∃ (E₀ : C) (ϵ₀ : E₀ ⟶ F), IsCoAngular P ϵ₀ ∧ IsPreStep P ϵ₀ ∧ IsPrimaryElt (P.Div ϵ₀) ∧
    ∀ (E' : C) (ϵ' : E' ⟶ F), IsCoAngular P ϵ' → IsPreStep P ϵ' → IsPrimaryElt (P.Div ϵ') →
      ¬ SamePrimeCat P ϵ' ϵ₀ →
        ((∃ ζ : E ⟶ E', IsPreStep P ζ ∧ ϵ = ζ ≫ ϵ') ↔
         (∃ θ : Do ⟶ E', IsPreStep P θ ∧ δ ≫ ϵ = θ ≫ ϵ'))

/-- ★★どの素点も `𝒪^▷` の元が定める primary pre-step で実現される。 -/
theorem exists_otri_realizing
    (hdivS : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
      ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
    {A : C} (p : Prime (Φ.val (P.toElem.obj A).base)) :
    ∃ (u : End A) (hu : u ∈ OTri P A) (h : IsPreStep P ((u : A ⟶ A)))
      (hp : IsPrimaryElt (preStepVal P ((u : A ⟶ A)) h)),
      toPrime _ (preStepVal P ((u : A ⟶ A)) h) hp = p := by
  obtain ⟨⟨a, ha⟩, hx⟩ := Quotient.exists_rep p
  obtain ⟨u, hu⟩ := hdivS A a
  have hval : preStepVal P (((u : End A)) : A ⟶ A) (isPreStep_of_otri _ u.2) = a := by
    rw [preStepVal_of_otri _ u.2]; exact hu
  subst hval
  exact ⟨(u : End A), u.2, isPreStep_of_otri _ u.2, ha, hx⟩

def Prop41ivRhsCat.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 77,
    item := "Proposition 4.1, (iv) — 右辺の圏論版",
    sectionId := "frdi-prop-4-1" }

/-- ★★★★★★**[FrdI] Proposition 4.1, (iv) の圏論版** ——
右辺から `Prime` と `SuppElt` が消える。 -/
theorem prop_4_1_iv_cat (G : Frobenioid P) (hiso : ∀ X : C, IsIsotropic P X)
    (hdivS : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
      ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
    {Do E F : C} (δ : Do ⟶ E) (ϵ : E ⟶ F)
    (hsδ : IsPreStep P δ) (hsϵ : IsPreStep P ϵ) (hδ0 : P.Div δ ≠ 0)
    {ι₀ : Prime (Φ.val (P.toElem.obj F).base) → Pf (Φ.val (P.toElem.obj F).base) → ℝ≥0}
    (H : IsPerfFactorialWith (Φ.val (P.toElem.obj F).base) ι₀)
    (hperf : IsPerfectMonoid (Φ.val (P.toElem.obj F).base))
    (hdivM : IsDivisorial (Φ.val (P.toElem.obj F).base)) :
    IsPrimaryElt (P.Div δ) ↔ Prop41ivRhsCat P δ ϵ := by
  rw [prop_4_1_iv G hiso δ ϵ hsδ hsϵ hδ0 H hperf]
  constructor
  · rintro ⟨p, hp⟩
    obtain ⟨u, hu, hus, hup, hutp⟩ := exists_otri_realizing hdivS p
    refine ⟨F, ((u : End F) : F ⟶ F), prop_1_4_i P _ (fun Y _ => hiso Y), hus,
      (isPrimaryElt_preStepVal_iff _ hus).mp hup, ?_⟩
    intro E' ϵ' hcϵ' hsϵ' hpϵ' hns
    refine hp E' ϵ' hsϵ' ((isPrimaryElt_preStepVal_iff ϵ' hsϵ').mpr hpϵ') ?_
    intro hsupp
    refine hns ?_
    refine (toPrime_preStepVal_eq_iff G hiso H hperf hdivM ϵ' hcϵ' hsϵ'
      ((isPrimaryElt_preStepVal_iff ϵ' hsϵ').mpr hpϵ')
      ((u : End F) : F ⟶ F) (prop_1_4_i P _ (fun Y _ => hiso Y)) hus hup).mp ?_
    have h1 : SuppElt ι₀ (preStepVal P ϵ' hsϵ')
        = {toPrime _ (preStepVal P ϵ' hsϵ')
            ((isPrimaryElt_preStepVal_iff ϵ' hsϵ').mpr hpϵ')} :=
      suppElt_eq_singleton_toPrime H hperf hdivM _
    rw [h1] at hsupp
    have h2 : toPrime _ (preStepVal P ϵ' hsϵ')
        ((isPrimaryElt_preStepVal_iff ϵ' hsϵ').mpr hpϵ') = p := Set.singleton_injective hsupp
    rw [h2, ← hutp]
  · rintro ⟨E₀, ϵ₀, hc₀, hs₀, hp₀, hcond⟩
    have hpv₀ : IsPrimaryElt (preStepVal P ϵ₀ hs₀) :=
      (isPrimaryElt_preStepVal_iff ϵ₀ hs₀).mpr hp₀
    refine ⟨toPrime _ (preStepVal P ϵ₀ hs₀) hpv₀, ?_⟩
    intro E' ϵ' hsϵ' hpv' hne
    have hcϵ' : IsCoAngular P ϵ' := prop_1_4_i P _ (fun Y _ => hiso Y)
    refine hcond E' ϵ' hcϵ' hsϵ' ((isPrimaryElt_preStepVal_iff ϵ' hsϵ').mp hpv') ?_
    intro hsame
    refine hne ?_
    have h3 : toPrime _ (preStepVal P ϵ' hsϵ') hpv'
        = toPrime _ (preStepVal P ϵ₀ hs₀) hpv₀ :=
      (toPrime_preStepVal_eq_iff G hiso H hperf hdivM ϵ' hcϵ' hsϵ' hpv'
        ϵ₀ hc₀ hs₀ hpv₀).mpr hsame
    rw [suppElt_eq_singleton_toPrime H hperf hdivM hpv', h3]

def prop_4_1_iv_cat.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 77,
    item := "Proposition 4.1, (iv) — 右辺の圏論版",
    sectionId := "frdi-prop-4-1" }

/-! ## ★同型の前置に対する不変性 -/

/-- ★★同型を前置しても `Div` の primary 性は変わらない。 -/
theorem isPrimaryElt_div_iso_comp {X Y Z : C} (e : X ⟶ Y) [IsIso e] (f : Y ⟶ Z)
    (hf : IsPreStep P f) :
    IsPrimaryElt (P.Div (e ≫ f)) ↔ IsPrimaryElt (P.Div f) := by
  have hd : P.Div (e ≫ f) = Φ.map (P.Base e) (P.Div f) := by
    rw [P.Div_comp, show P.Div e = 0 from isIsometric_of_isIso P e, smul_zero, add_zero]
  rw [hd]
  haveI : IsIso (P.Base e) := by
    rw [show P.Base e = P.Base e from rfl]
    exact ⟨⟨P.Base (inv e), by rw [← P.Base_comp, IsIso.hom_inv_id, P.Base_id],
      by rw [← P.Base_comp, IsIso.inv_hom_id, P.Base_id]⟩⟩
  constructor
  · intro hp
    have hbij : Function.Bijective (Φ.map (@inv _ _ _ _ (P.Base e) inferInstance)) :=
      Φ.map_bijective_of_iso (@asIso _ _ _ _ (P.Base e) inferInstance).symm
    have h2 := isPrimaryElt_of_bijective _ hbij hp
    rwa [MonoidOn.map_map_inv Φ (P.Base e) inferInstance (P.Div f)] at h2
  · exact isPrimaryElt_of_bijective (Φ.map (P.Base e))
      (Φ.map_bijective_of_iso (@asIso _ _ _ _ (P.Base e) inferInstance))

/-- ★★同型を前置しても `SamePrimeCat` は変わらない。 -/
theorem samePrimeCat_iso_comp {A B B' X : C} (e : X ⟶ B) [IsIso e] {ϵ : B ⟶ A} {ϵ' : B' ⟶ A} :
    SamePrimeCat P (e ≫ ϵ) ϵ' ↔ SamePrimeCat P ϵ ϵ' := by
  constructor
  · rintro ⟨Z, ζ, hζc, hζst, ⟨ψ, hψs, hfac⟩, hr⟩
    refine ⟨Z, ζ, hζc, hζst, ⟨inv e ≫ ψ, IsPreStep.comp P (isPreStep_of_isIso P (inv e)) hψs,
      ?_⟩, hr⟩
    rw [Category.assoc, ← hfac, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
  · rintro ⟨Z, ζ, hζc, hζst, ⟨ψ, hψs, hfac⟩, hr⟩
    refine ⟨Z, ζ, hζc, hζst, ⟨e ≫ ψ, IsPreStep.comp P (isPreStep_of_isIso P e) hψs, ?_⟩, hr⟩
    rw [Category.assoc, ← hfac]

/-! ## ★「pre-step を通る分解」の移送 -/

/-- ★`f` が `g` を通って pre-step で分解する。 -/
def FactorsPre (P : PreFrobenioid C Φ) {X Y Z : C} (f : X ⟶ Z) (g : Y ⟶ Z) : Prop :=
  ∃ ζ : X ⟶ Y, IsPreStep P ζ ∧ f = ζ ≫ g

/-- ★同型を前置しても変わらない。 -/
theorem factorsPre_iso_comp {X Y Y' Z : C} (e : Y ⟶ Y') [IsIso e] (f : X ⟶ Z) (g : Y' ⟶ Z) :
    FactorsPre P f (e ≫ g) ↔ FactorsPre P f g := by
  constructor
  · rintro ⟨ζ, hζ, hfac⟩
    exact ⟨ζ ≫ e, IsPreStep.comp P hζ (isPreStep_of_isIso P e), by rw [hfac, Category.assoc]⟩
  · rintro ⟨ζ, hζ, hfac⟩
    refine ⟨ζ ≫ inv e, IsPreStep.comp P hζ (isPreStep_of_isIso P (inv e)), ?_⟩
    rw [hfac, Category.assoc, ← Category.assoc (inv e), IsIso.inv_hom_id, Category.id_comp]

end Cat41iv

section Transport41iv

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}

/-- ★★「pre-step を通る分解」は圏同値で行き来する。 -/
theorem factorsPre_map_iff (Ψ : C ≌ C₂)
    (hPS : ∀ {X Y : C} (f : X ⟶ Y), IsPreStep P f → IsPreStep P₂ (Ψ.functor.map f))
    (hPS' : ∀ {X Y : C} (f : X ⟶ Y), IsPreStep P₂ (Ψ.functor.map f) → IsPreStep P f)
    {X Y Z : C} (f : X ⟶ Z) (g : Y ⟶ Z) :
    FactorsPre P f g ↔ FactorsPre P₂ (Ψ.functor.map f) (Ψ.functor.map g) := by
  constructor
  · rintro ⟨ζ, hζ, hfac⟩
    exact ⟨Ψ.functor.map ζ, hPS ζ hζ, by rw [hfac, Ψ.functor.map_comp]⟩
  · rintro ⟨ζ', hζ', hfac⟩
    obtain ⟨ζ, hζeq⟩ := Ψ.functor.map_surjective ζ'
    refine ⟨ζ, hPS' ζ (by rw [hζeq]; exact hζ'), ?_⟩
    refine Ψ.functor.map_injective ?_
    rw [Ψ.functor.map_comp, hζeq, hfac]

end Transport41iv

section Cat41ivEnd

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

end Cat41ivEnd

section MainTransport

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}

/-- ★★★★★★**`Proposition 4.1, (iv)` の右辺は `Ψ` で移る**。 -/
theorem prop41ivRhsCat_map (Ψ : C ≌ C₂)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hCA : ∀ {X Y : C} (f : X ⟶ Y), IsCoAngular P f → IsCoAngular P₂ (Ψ.functor.map f))
    (hPS : ∀ {X Y : C} (f : X ⟶ Y), IsPreStep P f → IsPreStep P₂ (Ψ.functor.map f))
    (hPS' : ∀ {X Y : C} (f : X ⟶ Y), IsPreStep P₂ (Ψ.functor.map f) → IsPreStep P f)
    (hST : ∀ {X Y : C} (f : X ⟶ Y), IsStep P f → IsStep P₂ (Ψ.functor.map f))
    {Do E F : C} {δ : Do ⟶ E} {ϵ : E ⟶ F}
    (hFwd : ∀ (E' : C) (ϵ' : E' ⟶ F), IsPreStep P ϵ' →
      IsPrimaryElt (P.Div ϵ') → IsPrimaryElt (P₂.Div (Ψ.functor.map ϵ')))
    (hBwd : ∀ (E' : C) (ϵ' : E' ⟶ F), IsPreStep P ϵ' →
      IsPrimaryElt (P₂.Div (Ψ.functor.map ϵ')) → IsPrimaryElt (P.Div ϵ'))
    (h : Prop41ivRhsCat P δ ϵ) :
    Prop41ivRhsCat P₂ (Ψ.functor.map δ) (Ψ.functor.map ϵ) := by
  obtain ⟨E₀, ϵ₀, hc₀, hs₀, hp₀, hcond⟩ := h
  refine ⟨Ψ.functor.obj E₀, Ψ.functor.map ϵ₀, hCA ϵ₀ hc₀, hPS ϵ₀ hs₀,
    hFwd E₀ ϵ₀ hs₀ hp₀, ?_⟩
  intro E' ϵ' hcϵ' hsϵ' hpϵ' hns
  obtain ⟨E'₀, ⟨e⟩⟩ := Functor.EssSurj.mem_essImage (F := Ψ.functor) E'
  obtain ⟨ϵ'₀, hϵ'₀⟩ := Ψ.functor.map_surjective (e.hom ≫ ϵ')
  have hsϵ'₀ : IsPreStep P ϵ'₀ := by
    refine hPS' ϵ'₀ ?_
    rw [hϵ'₀]
    exact IsPreStep.comp P₂ (isPreStep_of_isIso P₂ e.hom) hsϵ'
  have hpϵ'₀ : IsPrimaryElt (P.Div ϵ'₀) := by
    refine hBwd E'₀ ϵ'₀ hsϵ'₀ ?_
    rw [hϵ'₀]
    exact (isPrimaryElt_div_iso_comp e.hom ϵ' hsϵ').mpr hpϵ'
  have hns₀ : ¬ SamePrimeCat P ϵ'₀ ϵ₀ := by
    intro hc
    refine hns ?_
    have h1 := samePrimeCat_map Ψ.functor hCA hPS hST hc
    rw [hϵ'₀] at h1
    exact (samePrimeCat_iso_comp e.hom).mp h1
  have hiff := hcond E'₀ ϵ'₀ (prop_1_4_i P _ (fun Y _ => hiso Y)) hsϵ'₀ hpϵ'₀ hns₀
  have hL : FactorsPre P₂ (Ψ.functor.map ϵ) ϵ' ↔ FactorsPre P ϵ ϵ'₀ := by
    rw [factorsPre_map_iff Ψ hPS hPS' ϵ ϵ'₀, hϵ'₀]
    exact (factorsPre_iso_comp e.hom (Ψ.functor.map ϵ) ϵ').symm
  have hR : FactorsPre P₂ (Ψ.functor.map δ ≫ Ψ.functor.map ϵ) ϵ'
      ↔ FactorsPre P (δ ≫ ϵ) ϵ'₀ := by
    rw [factorsPre_map_iff Ψ hPS hPS' (δ ≫ ϵ) ϵ'₀, hϵ'₀, Ψ.functor.map_comp]
    exact (factorsPre_iso_comp e.hom (Ψ.functor.map δ ≫ Ψ.functor.map ϵ) ϵ').symm
  exact hL.trans (hiff.trans hR.symm)

def prop41ivRhsCat_map.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 79,
    item := "Theorem 4.2, (i) — Proposition 4.1, (iv) の右辺は Ψ で移る",
    sectionId := "frdi-thm-4-2" }

/-- ★★★★★★**primary 性の「一段の拡張」** ——
`F` へ入る pre-step について primary 性が両向きに移るなら、
`F` へ pre-step `ϵ : E ⟶ F` を持つ `E` へ入る pre-step についても移る。

★★これが原文 p.79 の「`A₁` → `B₁` → `C₁`」の各段である。 -/
theorem isPrimaryElt_div_extend (Ψ : C ≌ C₂)
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
    (hδ0 : P.Div δ ≠ 0) (hδ0₂ : P₂.Div (Ψ.functor.map δ) ≠ 0)
    {ι₀ : Prime (Φ.val (P.toElem.obj F).base) → Pf (Φ.val (P.toElem.obj F).base) → ℝ≥0}
    (H : IsPerfFactorialWith (Φ.val (P.toElem.obj F).base) ι₀)
    (hperf : IsPerfectMonoid (Φ.val (P.toElem.obj F).base))
    (hdivM : IsDivisorial (Φ.val (P.toElem.obj F).base))
    {ι₂ : Prime (Φ₂.val (P₂.toElem.obj (Ψ.functor.obj F)).base) →
      Pf (Φ₂.val (P₂.toElem.obj (Ψ.functor.obj F)).base) → ℝ≥0}
    (H₂ : IsPerfFactorialWith (Φ₂.val (P₂.toElem.obj (Ψ.functor.obj F)).base) ι₂)
    (hperf₂ : IsPerfectMonoid (Φ₂.val (P₂.toElem.obj (Ψ.functor.obj F)).base))
    (hdivM₂ : IsDivisorial (Φ₂.val (P₂.toElem.obj (Ψ.functor.obj F)).base))
    (hFwd : ∀ (E' : C) (ϵ' : E' ⟶ F), IsPreStep P ϵ' →
      IsPrimaryElt (P.Div ϵ') → IsPrimaryElt (P₂.Div (Ψ.functor.map ϵ')))
    (hBwd : ∀ (E' : C) (ϵ' : E' ⟶ F), IsPreStep P ϵ' →
      IsPrimaryElt (P₂.Div (Ψ.functor.map ϵ')) → IsPrimaryElt (P.Div ϵ'))
    (h : IsPrimaryElt (P.Div δ)) :
    IsPrimaryElt (P₂.Div (Ψ.functor.map δ)) :=
  (prop_4_1_iv_cat G₂ hiso₂ hdivS₂ (Ψ.functor.map δ) (Ψ.functor.map ϵ)
      (hPS δ hsδ) (hPS ϵ hsϵ) hδ0₂ H₂ hperf₂ hdivM₂).mpr
    (prop41ivRhsCat_map Ψ hiso hCA hPS hPS' hST hFwd hBwd
      ((prop_4_1_iv_cat G hiso hdivS δ ϵ hsδ hsϵ hδ0 H hperf hdivM).mp h))

def isPrimaryElt_div_extend.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 79,
    item := "Theorem 4.2, (i) — primary 性の一段の拡張",
    sectionId := "frdi-thm-4-2" }

end MainTransport

end ABC3.Found.FrdI
