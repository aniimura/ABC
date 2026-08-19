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

end Cat41iv

end ABC3.Found.FrdI
